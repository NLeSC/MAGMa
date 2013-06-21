"""Module for job submission and job data alteration/retrieval
"""
import uuid
import os
import csv
import StringIO
import json
import shutil
from sqlalchemy import create_engine, and_
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import sessionmaker, aliased, scoped_session
from sqlalchemy.sql import func
from sqlalchemy.sql.expression import desc, asc, null
from sqlalchemy.orm.exc import NoResultFound
import transaction
import requests
from magmaweb.models import Base, Metabolite, Scan, Peak, Fragment, Run
import magmaweb.user
from .jobquery import JobQuery


class ScanRequiredError(Exception):
    """Raised when a scan identifier is required, but non is supplied"""
    pass


class ScanNotFound(Exception):
    """Raised when a scan identifier is not found"""
    pass


class FragmentNotFound(Exception):
    """Raised when a fragment is not found"""
    pass


class JobIdException(Exception):
    def __init__(self, message, jobid):
        Exception.__init__(self, message)
        self.jobid = jobid


class JobException(Exception):
    def __init__(self,
                 job,
                 message='Calculation failed for an unknown reason',
                 ):
        Exception.__init__(self)
        self.message = message
        self.job = job


class JobNotFound(JobIdException):
    """Raised when a job with a identifier is not found"""


class JobSubmissionError(IOError):
    """Raised when a job fails to be submitted"""


class JobIncomplete(JobException):
    """Raised when a complete job is required but job isnt"""


class JobError(JobException):
    """Job which failed during run"""


class MissingDataError(JobError):
    """Raised when job is missing data like molecules, peaks, or fragments"""


def make_job_factory(params):
    """Returns :class:`JobFactory` instance based on ``params`` dict

    All keys starting with 'jobfactory.' are used.
    From keys 'jobfactory.' is removed.

    Can be used to create job factory from a config file
    """
    prefix = 'jobfactory.'
    d = {k.replace(prefix, ''): v
         for k, v in params.iteritems()
         if k.startswith(prefix)}
    return JobFactory(**d)


class JobFactory(object):
    """Factory which can create jobs """
    def __init__(self,
                 root_dir,
                 init_script='',
                 tarball=None,
                 script_fn='script.sh',
                 db_fn='results.db',
                 submit_url='http://localhost:9998',
                 ):
        """
        root_dir
            Directory in which jobs are created, retrieved

        db_fn
            Sqlite db file name in job directory (default results.db)

        submit_url
            Url of job launcher daemon where job can be submitted
            (default http://localhost:9998)

        script_fn
            Script in job directory which must be run by job launcher daemon

        init_script
            String containing os commands to put magma in path,
            eg. activate virtualenv or unpack tarball

        tarball
            Local absolute location of tarball which contains application
            which init_script will unpack and run
        """
        self.root_dir = root_dir
        self.db_fn = db_fn
        self.submit_url = submit_url
        self.script_fn = script_fn
        self.tarball = tarball
        self.init_script = init_script

    def _makeJobSession(self, jobid):
        """ Create job db connection """
        engine = create_engine(self.id2url(jobid))
        try:
            engine.connect()
        except OperationalError:
            raise JobNotFound("Data of job not found", jobid)
        return scoped_session(sessionmaker(bind=engine))

    def fromId(self, jobid):
        """Finds job db in job root dir.
        Returns a :class:`Job` instance

        ``jobid``
            Job identifier as uuid

        Raises :class:`JobNotFound` exception when job is not found by jobid
        """
        meta = self._getJobMeta(jobid)
        session = self._makeJobSession(jobid)
        db = JobDb(session)
        jdir = self.id2jobdir(jobid)
        return Job(meta, jdir, db)

    def _makeJobDir(self, jobid):
        """Creates job directory

        Returns path of directory
        """
        jdir = self.id2jobdir(jobid)
        os.makedirs(jdir)
        return jdir

    def _getJobMeta(self, jobid):
        """Returns JobMeta belonging to job with `jobid` identifier.

        Raises :class:`JobNotFound` exception when job is not found by jobid
        """
        meta = magmaweb.user.JobMeta.by_id(jobid)
        if meta is None:
            raise JobNotFound("Job not found in database", jobid)
        return meta

    def _addJobMeta(self, jobmeta):
        """Adds JobMeta to database"""
        magmaweb.user.JobMeta.add(jobmeta)

    def fromDb(self, dbfile, owner):
        """A job directory is created and the dbfile is copied into it
        Returns a Job instance

        ``dbfile``
            The sqlite result db

        ``owner``
            User id of owner

        Returns a :class:`Job` instance
        """
        jobid = uuid.uuid4()

        # create job dir
        jdir = self._makeJobDir(jobid)

        #copy dbfile into jobdir
        self._copyFile(dbfile, jobid)

        # session for job db
        session = self._makeJobSession(jobid)
        db = JobDb(session)
        run = db.runInfo()

        #register job in user db
        jobmeta = magmaweb.user.JobMeta(jobid, owner,
                                        description=run.description,
                                        ms_filename=run.ms_filename)
        self._addJobMeta(jobmeta)

        return Job(jobmeta, jdir, db)

    def fromScratch(self, owner):
        """Job from scratch with empty results db

        ``owner``
            User id of owner

        Returns a :class:`Job` instance
        """
        jobid = uuid.uuid4()

        # create job dir
        jdir = self._makeJobDir(jobid)

        # session for job db
        session = self._makeJobSession(jobid)
        Base.metadata.create_all(session.connection())  # @UndefinedVariable
        db = JobDb(session)

        #register job in user db
        jobmeta = magmaweb.user.JobMeta(jobid, owner)
        self._addJobMeta(jobmeta)

        return Job(jobmeta, jdir, db)

    def _copyFile(self, src, jid):
        """Copy content of file object 'src' to
        result db of job with identifier 'jid'.
        """
        src.seek(0)  # start from beginning
        dst = open(self.id2db(jid), 'wb')
        shutil.copyfileobj(src, dst, 2 << 16)
        dst.close()

    def cloneJob(self, job, owner):
        """Returns new Job which has copy of 'job's database.

        ``job``
            :class:`Job` job to clone

        ``owner``
            User id of owner

        Returns a :class:`Job` instance
        """
        jobid = uuid.uuid4()

        # create job dir
        jdir = self._makeJobDir(jobid)

        #copy db of old job into new jobdir
        src = open(self.id2db(job.id))
        self._copyFile(src, jobid)
        src.close()

        # session for job db
        session = self._makeJobSession(jobid)
        db = JobDb(session)

        #register job in user db
        jobmeta = magmaweb.user.JobMeta(jobid, owner,
                                        description=job.description,
                                        ms_filename=job.ms_filename,
                                        parentjobid=job.id)
        self._addJobMeta(jobmeta)

        return Job(jobmeta, jdir, db)

    def submitJob2Launcher(self, body):
        """Submits job query to job launcher daemon

        body is a dict which is submitted as json.

        returns job url of job submitted to job launcher
        (not the same as job.id).
        """
        headers = {'Content-Type': 'application/json',
                   'Accept': 'application/json',
                   }
        response = requests.post(self.submit_url,
                                 data=json.dumps(body),
                                 headers=headers,
                                 )

        return response.headers['Location']

    def submitQuery(self, query, job):
        """Writes job query script to job dir
        and submits job query to job launcher.

        Changes the job state to u'INITIAL' or
        u'SUBMISSION_ERROR' if job submission fails.
        `job.launcher_url` gets filled with url returned by job launcher.

        query is a :class:`JobQuery` object

        job is a :class:`Job` object

        Raises JobSubmissionError if job submission fails.
        """

        # write job script into job dir
        script_fn = os.path.join(query.dir, self.script_fn)
        script = open(script_fn, 'w')
        script.write(self.init_script)
        script.write("\n")  # hard to add newline in ini file so add it here
        script.write(query.script.format(db=self.db_fn, magma='magma'))
        script.close()

        body = {
            'jobdir': query.dir + '/',
            'executable': '/bin/sh',
            'prestaged': [self.script_fn, self.db_fn],
            'poststaged': [self.db_fn],
            'stderr': 'stderr.txt',
            'stdout': 'stdout.txt',
            'arguments': [self.script_fn],
            'status_callback_url': query.status_callback_url
        }
        body['prestaged'].extend(query.prestaged)

        if (self.tarball is not None):
            body['prestaged'].append(self.tarball)

        job.state = u'PENDING'

        # Job launcher will try to report states of the submitted job
        # even before this request has been handled
        # The job launcher will not be able to find the job
        # causing the first state updates not to be saved
        # because commits happen at the end of this request
        # So force commit now
        transaction.commit()
        # due to commit job.meta will become deattached
        # making job.meta unusable
        # by merging it the job.meta is attached and valid again
        job.meta = magmaweb.user.DBSession().merge(job.meta)

        try:
            launcher_url = self.submitJob2Launcher(body)
            # store launcher url so job can be cancelled later
            job.launcher_url = launcher_url
        except requests.exceptions.RequestException:
            job.state = u'SUBMISSION_ERROR'
            raise JobSubmissionError()

    def id2jobdir(self, jid):
        """Returns job directory based on jid and root_dir """
        return os.path.join(self.root_dir, str(jid))

    def id2db(self, jid):
        """Returns sqlite db of job with jid"""
        return os.path.join(self.id2jobdir(jid), self.db_fn)

    def id2url(self, jid):
        """Returns sqlalchemy url of sqlite db of job with jid """
        # 3rd / is for username:pw@host which sqlite does not need
        return 'sqlite:///' + self.id2db(jid)

    def cancel(self, job):
        """Cancel in-complete :class:`Job` on the job server
        """
        url = job.launcher_url
        return requests.delete(url)


class Job(object):
    """Job contains results database of Magma calculation run"""

    def __init__(self, meta, directory, db=None):
        """
        meta
            :class:`magmaweb.user.JobMeta` instance.
            For owner, parent etc.

        directory
            Directory where input and output files reside

        db
            :class:`JobDb` instance
        """
        self.dir = directory
        self.meta = meta

        if db is not None:
            self.db = db

    @property
    def id(self):
        """Identifier of the job, is a :class:`uuid.UUID`"""
        return self.meta.jobid

    @property
    def __name__(self):
        """Same as id and as lookup key
        using :class:`magmaweb.user.JobIdFactory`
        """
        return str(self.id)

    def jobquery(self, status_callback_url, restricted):
        """Returns :class:`JobQuery`

        'status_callback_url' is the url to PUT status of job to.
        'restricted' whether to build queries in restricted mode.
        """
        return JobQuery(self.dir,
                        status_callback_url=status_callback_url,
                        restricted=restricted,
                        )

    def get_description(self):
        """Description string of job"""
        return self.meta.description

    def set_description(self, description):
        """Sets description in JobMeta and JobDb.runInfo() if present"""
        self.meta.description = description
        run = self.db.runInfo()
        if run is not None:
            run.description = description

    description = property(get_description, set_description)

    def get_parent(self):
        """Identifier of parent job"""
        return self.meta.parentjobid

    def set_parent(self, parent):
        """Set parent of job"""
        self.meta.parentjobid = parent

    parent = property(get_parent, set_parent)

    def get_owner(self):
        """User id which owns this job"""
        return self.meta.owner

    def set_owner(self, userid):
        """Set owner of job"""
        self.meta.owner = userid

    owner = property(get_owner, set_owner)

    def get_state(self):
        """Returns state of job

        See ibis org.gridlab.gat.resources.Job.JobState for possible states.
        """
        return self.meta.state

    def set_state(self, newstate):
        """Set state of job"""
        self.meta.state = newstate

    state = property(get_state, set_state)

    def stderr(self):
        """Returns stderr text file or empty file if stderr does not exist"""
        try:
            return open(os.path.join(self.dir, 'stderr.txt'), 'rb')
        except IOError:
            return StringIO.StringIO()

    @property
    def created_at(self):
        """Datetime when job was created"""
        return self.meta.created_at

    def get_ms_filename(self):
        """Filename of MS datafile"""
        return self.meta.ms_filename

    def set_ms_filename(self, ms_filename):
        """Sets Filename of MS datafile in
        JobMeta and JobDb.runInfo(0 if present
        """
        self.meta.ms_filename = ms_filename
        run = self.db.runInfo()
        if run is not None:
            run.ms_filename = ms_filename

    ms_filename = property(get_ms_filename, set_ms_filename)

    def get_is_public(self):
        """Whether job is public or not"""
        return self.meta.is_public

    def set_is_public(self, is_public):
        """Set whether job is public or not

        True is public and False is private
        """
        self.meta.is_public = is_public

    is_public = property(get_is_public, set_is_public)

    def get_launcher_url(self):
        """Url of job on joblauncher"""
        return self.meta.launcher_url

    def set_launcher_url(self, launcher_url):
        """Set url of job on joblauncher"""
        self.meta.launcher_url = launcher_url

    launcher_url = property(get_launcher_url, set_launcher_url)

    def delete(self):
        """Deletes job from user database and deletes job directory"""
        magmaweb.user.JobMeta.delete(self.meta)
        # disconnect from job results database before removing database file
        self.db.session.remove()
        shutil.rmtree(self.dir)

    def is_complete(self, mustBeFilled=False):
        """Checks if job is complete

        If mustBeFilled==True then checks
            if jobs contains molecules, mspectras and fragments.

        Returns true or raises JobError or JobIncomplete or MissingDataError
        """
        if self.state == 'STOPPED':
            if mustBeFilled:
                if not self.db.hasMolecules():
                    raise MissingDataError(self, 'No molecules found')
                elif not self.db.hasMspectras():
                    raise MissingDataError(self, 'No mass spectras found')
                elif not self.db.hasFragments():
                    raise MissingDataError(self, 'No fragments found')
            return True
        elif self.state == 'ERROR':
            raise JobError(self)
        else:
            raise JobIncomplete(self)


class JobDb(object):
    """Database of a job"""
    def __init__(self, session):
        """SQLAlchemy session which is connected to database of job"""
        self.session = session

    def maxMSLevel(self):
        """Returns the maximum nr of MS levels """
        return self.session.query(func.max(Scan.mslevel)).scalar() or 0

    def runInfo(self):
        """Returns last run info or None if there is no run info"""
        # cache run info to prevent 'database is locked' errors
        if not hasattr(self, '_runInfo'):
            q = self.session.query(Run).order_by(Run.runid.desc())
            self._runInfo = q.first()
        return self._runInfo

    def metabolitesTotalCount(self):
        """Returns unfiltered and not paged count of metabolites"""
        return self.session.query(Metabolite).count()

    def extjsgridfilter(self, q, column, afilter):
        """Query helper to convert a extjs grid filter dict
        to a sqlalchemy query filter

        """
        if (afilter['type'] == 'numeric'):
            if (afilter['comparison'] == 'eq'):
                efilter = q.filter(column == afilter['value'])
            if (afilter['comparison'] == 'gt'):
                efilter = q.filter(column > afilter['value'])
            if (afilter['comparison'] == 'lt'):
                efilter = q.filter(column < afilter['value'])
        elif (afilter['type'] == 'string'):
            efilter = q.filter(column.contains(afilter['value']))
        elif (afilter['type'] == 'list'):
            efilter = q.filter(column.in_(afilter['value']))
        elif (afilter['type'] == 'boolean'):
            efilter = q.filter(column == afilter['value'])
        elif (afilter['type'] == 'null'):
            if not afilter['value']:
                efilter = q.filter(column == null())  # IS NULL
            else:
                efilter = q.filter(column != null())  # IS NOT NULL
        return efilter

    def _metabolitesQuery2Rows(self, start, limit, q):
        mets = []
        for r in q[start: limit + start]:
            met = r.Metabolite
            row = {'metid': met.metid,
                   'mol': met.mol,
                   'level': met.level,
                   'probability': met.probability,
                   'reactionsequence': met.reactionsequence,
                   'smiles': met.smiles,
                   'molformula': met.molformula,
                   'isquery': met.isquery,
                   'origin': met.origin,
                   'nhits': met.nhits,
                   'mim': met.mim,
                   'logp': met.logp,
                   'assigned': r.assigned > 0,
                   'reference': met.reference,
                   }
            if ('score' in r.keys()):
                row['score'] = r.score
                row['deltappm'] = r.deltappm
            mets.append(row)

        return mets

    def _addSortingToMetabolitesQuery(self,
                                      sorts,
                                      scanid,
                                      q,
                                      fragal,
                                      assign_q):
        for col3 in sorts:
            if col3['property'] == 'assigned':
                col2 = assign_q.c.assigned
            elif (col3['property'] == 'score'):
                if (scanid is not None):
                    col2 = fragal.score
                else:
                    raise ScanRequiredError()
            elif (col3['property'] == 'deltappm'):
                if (scanid is not None):
                    col2 = fragal.deltappm
                else:
                    raise ScanRequiredError()
            else:
                cprop = col3['property']
                col2 = Metabolite.__dict__[cprop]  # @UndefinedVariable
            if (col3['direction'] == 'DESC'):
                q = q.order_by(desc(col2))
            elif (col3['direction'] == 'ASC'):
                q = q.order_by(asc(col2))

        return q

    def metabolites(self,
                    start=0, limit=10,
                    sorts=None, scanid=None,
                    filters=None):
        """Returns dict with total and rows attribute

        start
            Offset
        limit
            Maximum nr of metabolites to return
        scanid
            Only return metabolites that have hits in scan with this identifier
            Adds score and deltappm columns
        filters
            List of dicts which is generated by
            ExtJS component Ext.ux.grid.FiltersFeature
        sorts
            How to sort metabolites. List of dicts. Eg.

        .. code-block:: python

                [{"property":"probability","direction":"DESC"},
                 {"property":"metid","direction":"ASC"}]

        """
        sorts = sorts or [{"property": "probability", "direction": "DESC"},
                          {"property": "metid", "direction": "ASC"}]
        filters = filters or []
        q = self.session.query(Metabolite)

        # custom filters
        fragal = aliased(Fragment)
        if (scanid is not None):
            # TODO: add score column + order by score
            q = q.add_columns(fragal.score, fragal.deltappm)
            q = q.join(fragal.metabolite)
            q = q.filter(fragal.parentfragid == 0)
            q = q.filter(fragal.scanid == scanid)

        # add assigned column
        assigned = func.count('*').label('assigned')
        assign_q = self.session.query(Peak.assigned_metid, assigned)
        assign_q = assign_q.filter(Peak.assigned_metid != null())
        assign_q = assign_q.group_by(Peak.assigned_metid).subquery()
        q = q.add_columns(assign_q.c.assigned).\
            outerjoin(assign_q, Metabolite.metid == assign_q.c.assigned_metid)

        for afilter in filters:
            if afilter['field'] == 'assigned':
                col = assign_q.c.assigned
                afilter['type'] = 'null'
            elif (afilter['field'] == 'score'):
                if (scanid is not None):
                    col = fragal.score
                else:
                    raise ScanRequiredError()
            elif (afilter['field'] == 'deltappm'):
                if (scanid is not None):
                    col = fragal.deltappm
                else:
                    raise ScanRequiredError()
            else:
                # generic filters
                ffield = afilter['field']
                col = Metabolite.__dict__[ffield]  # @UndefinedVariable
            q = self.extjsgridfilter(q, col, afilter)

        total = q.count()

        q = self._addSortingToMetabolitesQuery(sorts, scanid,
                                               q, fragal, assign_q)

        mets = self._metabolitesQuery2Rows(start, limit, q)

        return {'total': total, 'rows': mets}

    def metabolites2csv(self, metabolites, cols=None):
        """Converts array of metabolites to csv file handler

        Params:
          `metabolites`
            array like metabolites()['rows']
          `cols`
            Which metabolite columns should be returned.
            A empty list selects all columns.

        Return
        :class:`StringIO.StringIO`
        """
        cols = cols or []
        csvstr = StringIO.StringIO()
        headers = [
            'origin', 'smiles', 'probability', 'reactionsequence',
            'nhits', 'molformula', 'mim', 'isquery', 'logp',
            'reference'
        ]
        if ('score' in metabolites[0].keys()):
            headers.append('score')

        if (len(cols) > 0):
            crow = [key for key in cols if key in headers]
            headers = crow

        csvwriter = csv.DictWriter(csvstr, headers, extrasaction='ignore')
        csvwriter.writeheader()
        for m in metabolites:
            m['reactionsequence'] = '|'.join(m['reactionsequence'])
            csvwriter.writerow(m)

        return csvstr

    def metabolites2sdf(self, metabolites, cols=None):
        """Converts array of metabolites to a sdf string

        Params:
          `metabolites`
            array like metabolites()['rows']
          `cols`
            Which metabolite columns should be returned.
            A empty list selects all columns.

        Return
        String
        """
        s = ''
        cols = cols or []
        props = ['origin', 'smiles', 'probability', 'reactionsequence',
                 'nhits', 'molformula', 'mim', 'logp',
                 'reference']
        if ('score' in metabolites[0].keys()):
            props.append('score')

        if (len(cols) > 0):
            crow = [key for key in cols if key in props]
            props = crow

        for m in metabolites:
            s += m['mol']
            m['reactionsequence'] = '\n'.join(m['reactionsequence'])
            for p in props:
                s += '> <{}>\n{}\n\n'.format(p, m[p])
            s += '$$$$' + "\n"

        return s

    def scansWithMetabolites(self, filters=None, metid=None):
        """Returns id and rt of lvl1 scans which have a fragment in it
        and for which the filters in params pass

        params:

        ``metid``
            Only return scans that have hits with metabolite with this id
        ``filters``
            List of filters which is generated
            by ExtJS component Ext.ux.grid.FiltersFeature
        """
        filters = filters or []
        fq = self.session.query(Fragment.scanid).\
            filter(Fragment.parentfragid == 0)
        if (metid is not None):
            fq = fq.filter(Fragment.metid == metid)

        for afilter in filters:
            has_no_hit_filter = (afilter['field'] == 'nhits'
                                 and afilter['comparison'] == 'gt'
                                 and afilter['value'] == 0)
            if has_no_hit_filter:
                continue
            if (afilter['field'] == 'score'):
                fq = self.extjsgridfilter(fq, Fragment.score, afilter)
            elif (afilter['field'] == 'deltappm'):
                fq = self.extjsgridfilter(fq, Fragment.deltappm, afilter)
            elif (afilter['field'] == 'assigned'):
                afilter['type'] = 'null'
                fq = fq.join(Peak, and_(Fragment.scanid == Peak.scanid,
                                        Fragment.mz == Peak.mz))
                fq = self.extjsgridfilter(fq, Peak.assigned_metid, afilter)
            else:
                fq = fq.join(Metabolite,
                             Fragment.metabolite)  # @UndefinedVariable
                ffield = afilter['field']
                fcol = Metabolite.__dict__[ffield]  # @UndefinedVariable
                fq = self.extjsgridfilter(fq,
                                          fcol,
                                          afilter
                                          )

        hits = []
        q = self.session.query(Scan.rt, Scan.scanid).filter_by(mslevel=1)
        for hit in q.filter(Scan.scanid.in_(fq)):
            hits.append({
                'id': hit.scanid,
                'rt': hit.rt
            })

        return hits

    def extractedIonChromatogram(self, metid):
        """Returns extracted ion chromatogram of metabolite with id metid """
        chromatogram = []
        mzqq = self.session.query(func.avg(Fragment.mz))
        mzqq = mzqq.filter(Fragment.metid == metid)
        mzqq = mzqq.filter(Fragment.parentfragid == 0)
        mzq = mzqq.scalar()
        precision = 1 + self.session.query(Run.mz_precision).scalar() / 1e6
        # fetch max intensity of peaks with mz = mzq+-mzoffset
        q = self.session.query(Scan.rt, func.max(Peak.intensity))
        q = q.outerjoin(Peak, and_(Peak.scanid == Scan.scanid,
                                   Peak.mz.between(mzq / precision,
                                                   mzq * precision)))
        q = q.filter(Scan.mslevel == 1)
        for (rt, intens) in q.group_by(Scan.rt).order_by(asc(Scan.rt)):
            chromatogram.append({
                'rt': rt,
                'intensity': intens or 0
            })
        return chromatogram

    def chromatogram(self):
        """Returns dict with `cutoff` key with ms_intensity_cutoff and
        `scans` key which is a list of all lvl1 scans.

        Each scan is a dict with:
            * id, scan identifier
            * rt, retention time
            * intensity
        """
        scans = []

        assigned_peaks = func.count('*').label('assigned_peaks')
        ap = self.session.query(Peak.scanid, assigned_peaks)
        ap = ap.filter(Peak.assigned_metid != null())
        ap = ap.group_by(Peak.scanid).subquery()

        q = self.session.query(Scan, ap.c.assigned_peaks)
        q = q.filter_by(mslevel=1)
        for scan, assigned_peaks in q.outerjoin(ap,
                                                Scan.scanid == ap.c.scanid):
            scans.append({
                'id': scan.scanid,
                'rt': scan.rt,
                'intensity': scan.basepeakintensity,
                'ap': assigned_peaks or 0
            })

        runInfo = self.runInfo()
        if (runInfo is not None):
            return {'scans': scans, 'cutoff': runInfo.ms_intensity_cutoff}
        else:
            return {'scans': scans, 'cutoff': None}

    def mspectra(self, scanid, mslevel=None):
        """Returns dict with peaks of a scan

        Also returns the cutoff applied to the scan
        and mslevel, precursor.id (parent scan id) and precursor.mz

        scanid
            Scan identifier of scan of which to return the mspectra

        mslevel
            Ms level on which the scan must be. Optional.
            If scanid not on mslevel raises ScanNotFound

        """
        scanq = self.session.query(Scan).filter(Scan.scanid == scanid)
        if (mslevel is not None):
            scanq = scanq.filter(Scan.mslevel == mslevel)

        try:
            scan = scanq.one()
        except NoResultFound:
            raise ScanNotFound()

        # lvl1 scans use absolute cutoff, lvl>1 use ratio of basepeak as cutoff
        if (scan.mslevel == 1):
            cutoff = self.session.query(Run.ms_intensity_cutoff).scalar()
        else:
            rel_cutoff = Scan.basepeakintensity * Run.msms_intensity_cutoff
            q = self.session.query(rel_cutoff)
            cutoff = q.filter(Scan.scanid == scanid).scalar()

        peaks = []
        for peak in self.session.query(Peak).filter_by(scanid=scanid):
            peaks.append({
                'mz': peak.mz,
                'intensity': peak.intensity,
                'assigned_metid': peak.assigned_metid,
            })

        precursor = {'id': scan.precursorscanid, 'mz': scan.precursormz}
        return {'peaks': peaks,
                'cutoff': cutoff,
                'mslevel': scan.mslevel,
                'precursor': precursor,
                }

    def _fragmentsQuery(self):
        return self.session.query(Fragment,
                                  Metabolite.mol,
                                  Scan.mslevel).join(Metabolite).join(Scan)

    def _fragment2json(self, row):
        (frag, mol, mslevel) = row
        f = {
            'fragid': frag.fragid,
            'scanid': frag.scanid,
            'metid': frag.metid,
            'score': frag.score,
            'mol': mol,
            'atoms': frag.atoms,
            'mz': frag.mz,
            'mass': frag.mass,
            'deltah': frag.deltah,
            'mslevel': mslevel,
            'deltappm': frag.deltappm,
            'formula': frag.formula,
        }
        if (len(frag.children) > 0):
            f['expanded'] = False
            f['leaf'] = False
        else:
            f['expanded'] = True
            f['leaf'] = True
        return f

    def fragments(self, scanid, metid, node):
        """Returns fragments of a metabolite on a scan.

        When node is not set then returns metabolites and its lvl2 fragments.

        When node is set then returns children fragments as list
        which have ``node`` as parent fragment.

        Can be used in a Extjs.data.TreeStore if jsonified.

        Parameters:

        * ``scanid``, Fragments on scan with this identifier
        * ``metid``, Fragments of metabolite with this identifier
        * ``node``, The fragment identifier to fetch children fragments for.
            Use 'root' for root node.

        Raises FragmentNotFound when no fragment is found
        with the given scanid/metid combination
        """

        # parent metabolite
        if (node == 'root'):
            structures = []
            pms = self._fragmentsQuery().filter(Fragment.scanid == scanid)
            pms = pms.filter(Fragment.metid == metid)
            pms = pms.filter(Fragment.parentfragid == 0)
            for row in pms:
                structure = self._fragment2json(row)

                qa = self.session.query(func.count('*')).\
                    filter(Peak.scanid == scanid).\
                    filter(Peak.assigned_metid == metid)
                structure['isAssigned'] = qa.scalar() > 0

                # load children
                structure['children'] = []
                pf = Fragment.parentfragid
                cq = self._fragmentsQuery().filter(pf == structure['fragid'])
                for frow in cq:
                    structure['expanded'] = True
                    structure['children'].append(self._fragment2json(frow))
                structures.append(structure)

            if (len(structures) == 0):
                raise FragmentNotFound()
            return {'children': structures, 'expanded': True}
        # fragments
        else:
            fragments = []
            fq = self._fragmentsQuery().filter(Fragment.parentfragid == node)
            for row in fq:
                fragments.append(self._fragment2json(row))
            return fragments

    def _peak(self, scanid, mz):
        mzoffset = 1e-6  # precision for comparing floating point values
        # fetch peak corresponding to given scanid and mz
        q = self.session.query(Peak).filter(Peak.scanid == scanid)
        q = q.filter(Peak.mz.between(float(mz) - mzoffset,
                                     float(mz) + mzoffset))
        return q.one()

    def assign_metabolite2peak(self, scanid, mz, metid):
        """Assign metabolites to peak"""
        peak = self._peak(scanid, mz)
        peak.assigned_metid = metid
        self.session.add(peak)
        self.session.commit()

    def unassign_metabolite2peak(self, scanid, mz):
        """Unassign any metabolite from peak"""
        peak = self._peak(scanid, mz)
        peak.assigned_metid = None
        self.session.add(peak)
        self.session.commit()

    def hasMolecules(self):
        """Does job database contain molecules"""
        return self.session.query(Metabolite).count() > 0

    def hasMspectras(self):
        """Does job database contain mass spectra peaks"""
        return self.session.query(Peak).count() > 0

    def hasFragments(self):
        """Does job database contain fragments"""
        return self.session.query(Fragment).count() > 0
