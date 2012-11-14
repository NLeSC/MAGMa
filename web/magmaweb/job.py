import uuid
import os
import csv
import StringIO
import urllib2
import json
import shutil
from cgi import FieldStorage
import colander
from sqlalchemy import create_engine, and_
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import sessionmaker, aliased, scoped_session
from sqlalchemy.sql import func
from sqlalchemy.sql.expression import desc, asc, null
from sqlalchemy.orm.exc import NoResultFound
from magmaweb.models import Base, Metabolite, Scan, Peak, Fragment, Run
import magmaweb.user


class ScanRequiredError(Exception):
    """Raised when a scan identifier is required, but non is supplied"""
    pass


class ScanNotFound(Exception):
    """Raised when a scan identifier is not found"""
    pass


class FragmentNotFound(Exception):
    """Raised when a fragment is not found"""
    pass


class JobNotFound(Exception):
    """Raised when a job with a identifier is not found"""

    def __init__(self, message, jobid):
        Exception.__init__(self, message)
        self.jobid = jobid


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


class JobQuery(object):
    """Perform actions on a :class:`Job` by parsing post params,
    staging files and building magma cli commands
    """

    class File(object):
        """ Colander schema type for file upload field"""
        def serialize(self, node, appstruct):
            return appstruct

        def deserialize(self, node, cstruct):
            if cstruct is colander.null:
                return colander.null
            if cstruct == '':
                return colander.null
            if not isinstance(cstruct, FieldStorage):
                raise colander.Invalid(
                    node,
                    '{} is not a cgi.FieldStorage'.format(cstruct)
                )
            return cstruct

    def __init__(self, id, dir, script='', prestaged=None):
        """Contruct JobQuery

        Params:
        - id, job identifier
        - dir, job directory
        - script, script to run in job directory
        - prestaged, list of files to prestage

        """
        self.id = id
        self.dir = dir
        self.script = script
        self.prestaged = prestaged or []

    def __eq__(self, other):
        return (self.id == other.id and
                self.dir == other.dir and
                self.script == other.script and
                self.prestaged == other.prestaged)

    def __repr__(self):
        """Return a printable representation."""
        return "JobQuery({!r}, {!r}, {!r}, {!r})".format(self.id,
                                                         self.dir,
                                                         self.script,
                                                         self.prestaged)

    def escape(self, string):
        """ Replaces single quote with its html escape sequence"""
        return str(string).replace("'", '&#39;')

    def _addAnnotateSchema(self, schema):
        schema.add(colander.SchemaNode(colander.Float(),
                                       name='precursor_mz_precision'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       name='mz_precision'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       name='ms_intensity_cutoff'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       name='msms_intensity_cutoff'))
        schema.add(colander.SchemaNode(
                                       colander.Integer(),
                                       validator=colander.OneOf([-1, 1]),
                                       name='ionisation_mode'))
        schema.add(colander.SchemaNode(colander.Integer(),
                                       validator=colander.Range(min=0),
                                       name='max_broken_bonds'))
        schema.add(colander.SchemaNode(colander.Boolean(),
                                       default=False, missing=False,
                                       name='use_all_peaks'))
        schema.add(colander.SchemaNode(colander.Boolean(),
                                       default=False, missing=False,
                                       name='skip_fragmentation'))

    def _addMetabolizeSchema(self, schema):
        validator = colander.OneOf(['phase1', 'phase2'])
        metabolism_type = colander.SchemaNode(colander.String(),
                                              validator=validator)
        schema.add(colander.SchemaNode(colander.Sequence(),
                                       metabolism_type,
                                       name='metabolism_types'))
        schema.add(colander.SchemaNode(colander.Integer(),
                                       validator=colander.Range(0, 10),
                                       name='n_reaction_steps'))

    def add_structures(self, params, has_ms_data=False):
        """Configure job query to add_structures from params.

        ``params`` is a MultiDict from which the following keys are used:

        * structure_format, in which format is structures or structure_file
        * structures, string with structures
        * structures_file, file-like object with structures
        * metabolize, when key is set then
            :meth:`~magmaweb.job.JobQuery.metabolize` will be called.

        If ``has_ms_data`` is True then
            :meth:`~magmaweb.job.JobQuery.annotate` will be called.
        If both ``stuctures`` and ``structures_file`` is filled then
            ``structures`` is ignored.
        """
        def textarea_or_file(node, value):
            """Validator that either textarea or file upload is filled"""
            if not(not value['structures'] is colander.null or
                   not value['structures_file'] is colander.null):
                error = 'Either structures or structure_file must be set'
                exception = colander.Invalid(node)
                exception.add(colander.Invalid(structuresSchema, error))
                exception.add(colander.Invalid(structures_fileSchema, error))
                raise exception

        schema = colander.SchemaNode(colander.Mapping(),
                                     validator=textarea_or_file)
        schema.add(colander.SchemaNode(colander.String(),
                                       validator=colander.OneOf(['smiles',
                                                                 'sdf']),
                                       name='structure_format'
                                       ))
        filled = colander.Length(min=1)
        structuresSchema = colander.SchemaNode(colander.String(),
                                               validator=filled,
                                               missing=colander.null,
                                               name='structures')
        schema.add(structuresSchema)
        structures_fileSchema = colander.SchemaNode(self.File(),
                                                    missing=colander.null,
                                                    name='structures_file')
        schema.add(structures_fileSchema)
        schema.add(colander.SchemaNode(colander.Boolean(),
                                       default=False,
                                       missing=False,
                                       name='metabolize'))
        if has_ms_data:
            self._addAnnotateSchema(schema)
        if ('metabolize' in params):
            self._addMetabolizeSchema(schema)
        if hasattr(params, 'mixed'):
            # unflatten multidict
            params = params.mixed()
        if ('metabolism_types' in params and
                isinstance(params['metabolism_types'], basestring)):
            params['metabolism_types'] = [params['metabolism_types']]
        params = schema.deserialize(params)

        metsfile = file(os.path.join(self.dir, 'structures.dat'), 'w')
        if not params['structures_file'] is colander.null:
            sf = params['structures_file'].file
            sf.seek(0)
            while 1:
                data = sf.read(2 << 16)
                if not data:
                    break
                metsfile.write(data)
            sf.close()
        else:
            metsfile.write(params['structures'])
        metsfile.close()

        script = "{{magma}} add_structures -t '{structure_format}'"
        script += " structures.dat {{db}}"
        sf = self.escape(params['structure_format'])
        self.script += script.format(structure_format=sf)
        self.prestaged.append('structures.dat')

        if params['metabolize']:
            self.script += " |"
            self.metabolize(params, has_ms_data, True)
        elif (has_ms_data):
            # dont annotate when metabolize is also done
            # as metabolize will call annotate
            self.script += " |"
            self.annotate(params, True)
        else:
            self.script += "\n"

        return self

    def add_ms_data(self, params, has_metabolites=False):
        """Configure job query to add ms data from params.

        ``params`` is a dict from which the following keys are used:

        * ms_data_format
        * ms_data_file
        * max_ms_level
        * abs_peak_cutoff
        * rel_peak_cutoff

        If ``has_metabolites`` is True then
            :meth:`~magmaweb.job.JobQuery.annotate` will be called.

        """
        schema = colander.SchemaNode(colander.Mapping())
        schema.add(colander.SchemaNode(colander.String(),
                                       validator=colander.OneOf(['mzxml']),
                                       name='ms_data_format'
                                       ))
        schema.add(colander.SchemaNode(self.File(), name='ms_data_file'))
        schema.add(colander.SchemaNode(colander.Integer(),
                                       validator=colander.Range(min=0),
                                       name='max_ms_level'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       validator=colander.Range(min=0),
                                       name='abs_peak_cutoff'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       validator=colander.Range(0, 1),
                                       name='rel_peak_cutoff'))
        if has_metabolites:
            self._addAnnotateSchema(schema)
        params = schema.deserialize(params)

        msfile = file(os.path.join(self.dir, 'ms_data.dat'), 'w')
        msf = params['ms_data_file'].file
        msf.seek(0)
        while 1:
            data = msf.read(2 << 16)
            if not data:
                break
            msfile.write(data)
        msf.close()
        msfile.close()

        script = "{{magma}} read_ms_data --ms_data_format '{ms_data_format}' "
        script += "-l '{max_ms_level}' "
        script += "-a '{abs_peak_cutoff}' -r '{rel_peak_cutoff}' "
        script += "ms_data.dat {{db}}\n"
        script__substitution = {
            'ms_data_format': self.escape(params['ms_data_format']),
            'max_ms_level': self.escape(params['max_ms_level']),
            'abs_peak_cutoff': self.escape(params['abs_peak_cutoff']),
            'rel_peak_cutoff': self.escape(params['rel_peak_cutoff'])
        }
        self.script += script.format(**script__substitution)

        self.prestaged.append('ms_data.dat')

        if (has_metabolites):
            self.annotate(params)

        return self

    def metabolize(self, params, has_ms_data=False, from_subset=False):
        """Configure job query to metabolize all structures.

        ``params`` is a MultiDict from which the following keys are used:

        * n_reaction_steps
        * metabolism_types, comma seperated string with metabolism types

        If ``has_ms_data`` is True then
            :meth:`~magmaweb.job.JobQuery.annotate` will be called.
        If ``from_subset`` is True then metids are read from stdin
        """
        schema = colander.SchemaNode(colander.Mapping())
        self._addMetabolizeSchema(schema)
        if has_ms_data:
            self._addAnnotateSchema(schema)
        if hasattr(params, 'mixed'):
            # unflatten multidict
            params = params.mixed()
        if ('metabolism_types' in params and
                isinstance(params['metabolism_types'], basestring)):
            params['metabolism_types'] = [params['metabolism_types']]
        params = schema.deserialize(params)

        script = "{{magma}} metabolize -s '{n_reaction_steps}'"
        script += " -m '{metabolism_types}'"
        metabolism_types = ','.join(params['metabolism_types'])
        script_substitution = {
            'n_reaction_steps': self.escape(params['n_reaction_steps']),
            'metabolism_types': self.escape(metabolism_types)
        }
        self.script += script.format(**script_substitution)

        if from_subset:
            self.script += " -j -"

        if (has_ms_data):
            self.script += " {db} |"
            self.annotate(params, True)
        else:
            self.script += " {db}\n"

        return self

    def metabolize_one(self, params, has_ms_data=False):
        """Configure job query to metabolize one structure.

        ``params`` is a MultiDict from which the following keys are used:

        * metid, metabolite identifier to metabolize
        * n_reaction_steps
        * metabolism_types, comma seperated string with metabolism types

        If ``has_ms_data`` is True then
            :meth:`~magmaweb.job.JobQuery.annotate` will be called.
        """
        schema = colander.SchemaNode(colander.Mapping())
        schema.add(colander.SchemaNode(colander.Integer(),
                                       validator=colander.Range(min=0),
                                       name='metid'))
        self._addMetabolizeSchema(schema)
        if has_ms_data:
            self._addAnnotateSchema(schema)
        if hasattr(params, 'mixed'):
            # unflatten multidict
            params = params.mixed()
        if ('metabolism_types' in params and
                isinstance(params['metabolism_types'], basestring)):
            params['metabolism_types'] = [params['metabolism_types']]
        params = schema.deserialize(params)

        script = "echo '{metid}' | {{magma}} metabolize -j - "
        script += "-s '{n_reaction_steps}' -m '{metabolism_types}' {{db}}"
        metabolism_types = ','.join(params['metabolism_types'])
        script_substitutions = {
             'metid': self.escape(params['metid']),
             'n_reaction_steps': self.escape(params['n_reaction_steps']),
             'metabolism_types': self.escape(metabolism_types)
        }
        self.script += script.format(**script_substitutions)

        if (has_ms_data):
            self.script += " |"
            self.annotate(params, True)
        else:
            self.script += "\n"

        return self

    def annotate(self, params, from_subset=False):
        """Configure job query to annotate.

        ``params`` is a dict from which the following keys are used:

        * precursor_mz_precision
        * mz_precision
        * ms_intensity_cutoff
        * msms_intensity_cutoff
        * ionisation_mode
        * max_broken_bonds
        * use_all_peaks, when key is set then all peaks are used
        * skip_fragmentation, when key is set then
            the no fragmentation of structures is performed.

        If ``from_subset`` is True then metids are read from stdin
        """
        schema = colander.SchemaNode(colander.Mapping())
        self._addAnnotateSchema(schema)
        params = schema.deserialize(params)

        script = "{{magma}} annotate -p '{mz_precision}'"
        script += " -c '{ms_intensity_cutoff}' -d '{msms_intensity_cutoff}' "
        script += "-i '{ionisation_mode}' -b '{max_broken_bonds}'"
        script += " --precursor_mz_precision '{precursor_mz_precision}' "
        pmzp = params['precursor_mz_precision']
        ms_ic = params['ms_intensity_cutoff']
        msms_ic = params['msms_intensity_cutoff']
        script_substitutions = {
               'precursor_mz_precision': self.escape(pmzp),
               'mz_precision': self.escape(params['mz_precision']),
               'ms_intensity_cutoff': self.escape(ms_ic),
               'msms_intensity_cutoff': self.escape(msms_ic),
               'ionisation_mode': self.escape(params['ionisation_mode']),
               'max_broken_bonds': self.escape(params['max_broken_bonds'])
        }
        script = script.format(**script_substitutions)

        if (params['use_all_peaks']):
            script += '-u '

        if (params['skip_fragmentation']):
            script += '-f '

        if from_subset:
            script += '-j - '

        script += "{db}\n"
        self.script += script

        return self

    def allinone(self, params):
        """Configure job query to do all sub commands in one go.

        params is a MultiDict

        See
        :meth:`~magmaweb.job.JobQuery.add_ms_data`,
        :meth:`~magmaweb.job.JobQuery.add_structures`,
        :meth:`~magmaweb.job.JobQuery.metabolize`,
        :meth:`~magmaweb.job.JobQuery.annotate` for params.

        """
        allin = self.add_ms_data(params).add_structures(params)
        return allin.metabolize(params).annotate(params)

    @classmethod
    def defaults(cls):
        """Returns dictionary with default params"""
        return dict(
                    n_reaction_steps=2,
                    metabolism_types=['phase1', 'phase2'],
                    ionisation_mode=1,
                    skip_fragmentation=False,
                    ms_intensity_cutoff=1000000.0,
                    msms_intensity_cutoff=0.1,
                    mz_precision=0.001,
                    use_all_peaks=False,
                    abs_peak_cutoff=1000,
                    rel_peak_cutoff=0.01,
                    max_ms_level=10,
                    precursor_mz_precision=0.005,
                    max_broken_bonds=4
                    )


class JobFactory(object):
    """Factory which can create jobs """
    def __init__(
                  self, root_dir,
                  init_script='', tarball=None, script_fn='script.sh',
                  db_fn='results.db', state_fn='job.state',
                  submit_url='http://localhost:9998', time_max=30):
        """
        root_dir
            Directory in which jobs are created, retrieved

        db_fn
            Sqlite db file name in job directory (default results.db)

        submit_url
            Url of job manager daemon where job can be submitted
            (default http://localhost:9998)

        script_fn
            Script in job directory which must be run by job manager daemon

        init_script
            String containing os commands to put magma in path,
            eg. activate virtualenv or unpack tarball

        tarball
            Local absolute location of tarball which contains application
            which init_script will unpack and run

        state_fn
            Filename where job manager daemon writes job state
            (default job.state)

        time_max
            Maximum time in minutes a job can take (default 30)
        """
        self.root_dir = root_dir
        self.db_fn = db_fn
        self.submit_url = submit_url
        self.state_fn = state_fn
        self.script_fn = script_fn
        self.tarball = tarball
        self.time_max = time_max
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
        jdir = self.id2jobdir(jobid)
        os.makedirs(jdir)
        return jdir

    def _getJobMeta(self, jobid):
        meta = magmaweb.user.JobMeta.by_id(jobid)
        if meta is None:
            raise JobNotFound("Job not found in database", jobid)
        return meta

    def _addJobMeta(self, jobmeta):
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
        Base.metadata.create_all(session.connection())
        db = JobDb(session)

        #register job in user db
        jobmeta = magmaweb.user.JobMeta(jobid, owner)
        self._addJobMeta(jobmeta)

        return Job(jobmeta, jdir, db)

    def _copyFile(self, src, jid):
        """Copy content of file object 'src' to result db of job with identifier 'jid'."""
        src.seek(0)  # start from beginning
        dst = open(self.id2db(jid), 'wb')
        shutil.copyfileobj(src, dst, 2 << 16)
        dst.close()

    def cloneJob(self, job, owner):
        """Returns new Job which has copy of 'job's database.

        ``job``
            :class:Job job to clone

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

    def submitJob2Manager(self, body):
        """Submits job query to jobmanager daemon

        body is a dict which is submitted as json.

        returns job identifier of job submitted to jobmanager
        (not the same as job.id).
        """
        request = urllib2.Request(
                                  self.submit_url,
                                  json.dumps(body),
                                  {'Content-Type': 'application/json'})
        # log what is send to job manager
        import logging
        logger = logging.getLogger('magmaweb')
        logger.info(request.data)
        return urllib2.urlopen(request)

    def submitQuery(self, query):
        """Writes job script to job dir and submits job to job manager

        query is a :class:`JobQuery` object

        Returns job identifier
        """

        # write job script into job dir
        script = open(os.path.join(query.dir, self.script_fn), 'w')
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
               'time_max': self.time_max,
               'arguments': [self.script_fn]
               }
        body['prestaged'].extend(query.prestaged)

        if (self.tarball is not None):
            body['prestaged'].append(self.tarball)

        self.submitJob2Manager(body)

        return query.id

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


class Job(object):
    """Job contains results database of Magma calculation run"""

    def __init__(self, meta, dir, db=None):
        """
        meta
            :class:magmaweb.user.JobMeta instance.
            For owner, parent etc.

        dir
            Directory where input and output files reside

        db
            :class:JobDb instance
        """
        self.dir = dir
        self.meta = meta

        if db is not None:
            self.db = db

    @property
    def id(self):
        """Identifier of the job, is a :class:uuid.UUID"""
        return self.meta.jobid

    @property
    def __name__(self):
        """Same as id and as lookup key using :class:magmaweb.user.JobIdFactory"""
        return str(self.id)

    def jobquery(self):
        """Returns :class:`JobQuery` """
        return JobQuery(self.id, self.dir)

    @property
    def description(self):
        """Description string of job"""
        return self.meta.description

    @description.setter
    def description(self, description):
        """Sets description in JobMeta and JobDb.runInfo() if present"""
        self.meta.description = description
        run = self.db.runInfo()
        if run is not None:
            run.description = description

    @property
    def parent(self):
        """Identifier of parent job"""
        return self.meta.parentjobid

    @parent.setter
    def parent(self, parent):
        self.meta.parentjobid = parent

    @property
    def owner(self):
        """User id which owns this job"""
        return self.meta.owner

    @owner.setter
    def owner(self, userid):
        self.meta.owner = userid

    @property
    def state(self):
        """Returns state of job

        See ibis org.gridlab.gat.resources.Job.JobState for possible states.
        """
        return self.meta.state

    @state.setter
    def state(self, newstate):
        self.meta.state = newstate

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

    @property
    def ms_filename(self):
        """Filename of MS datafile"""
        return self.meta.ms_filename

    @ms_filename.setter
    def ms_filename(self, ms_filename):
        """Sets Filename of MS datafile in JobMeta and JobDb.runInfo(0 if present"""
        self.meta.ms_filename = ms_filename
        run = self.db.runInfo()
        if run is not None:
            run.ms_filename = ms_filename


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
        return self.session.query(Run).order_by(Run.runid.desc()).first()

    def metabolitesTotalCount(self):
        """Returns unfiltered and not paged count of metabolites"""
        return self.session.query(Metabolite).count()

    def extjsgridfilter(self, q, column, filter):
        """Query helper to convert a extjs grid filter dict
        to a sqlalchemy query filter

        """
        if (filter['type'] == 'numeric'):
            if (filter['comparison'] == 'eq'):
                return q.filter(column == filter['value'])
            if (filter['comparison'] == 'gt'):
                return q.filter(column > filter['value'])
            if (filter['comparison'] == 'lt'):
                return q.filter(column < filter['value'])
        elif (filter['type'] == 'string'):
            return q.filter(column.contains(filter['value']))
        elif (filter['type'] == 'list'):
            return q.filter(column.in_(filter['value']))
        elif (filter['type'] == 'boolean'):
            return q.filter(column == filter['value'])
        elif (filter['type'] == 'null'):
            if not filter['value']:
                return q.filter(column == null())  # IS NULL
            else:
                return q.filter(column != null())  # IS NOT NULL

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
        sorts = sorts or [{"property":"probability", "direction":"DESC"},
                          {"property":"metid", "direction":"ASC"}]
        filters = filters or []
        mets = []
        q = self.session.query(Metabolite)

        # custom filters
        fragal = aliased(Fragment)
        if (scanid is not None):
            # TODO: add score column + order by score
            q = q.add_columns(fragal.score, fragal.deltappm)
            q = q.join(fragal.metabolite).filter(
                fragal.parentfragid == 0).filter(fragal.scanid == scanid)

        # add assigned column
        assigned = func.count('*').label('assigned')
        stmt2 = self.session.query(Peak.assigned_metid, assigned)
        stmt2 = stmt2.filter(Peak.assigned_metid != null())
        stmt2 = stmt2.group_by(Peak.assigned_metid).subquery()
        q = q.add_columns(stmt2.c.assigned).\
            outerjoin(stmt2, Metabolite.metid == stmt2.c.assigned_metid)

        for filter in filters:
            if filter['field'] == 'assigned':
                col = stmt2.c.assigned
                filter['type'] = 'null'
            elif (filter['field'] == 'score'):
                if (scanid is not None):
                    col = fragal.score
                else:
                    raise ScanRequiredError()
            elif (filter['field'] == 'deltappm'):
                if (scanid is not None):
                    col = fragal.deltappm
                else:
                    raise ScanRequiredError()
            else:
                # generic filters
                col = Metabolite.__dict__[filter['field']]  #@UndefinedVariable @IgnorePep8
            q = self.extjsgridfilter(q, col, filter)

        total = q.count()

        for col in sorts:
            if col['property'] == 'assigned':
                col2 = stmt2.c.assigned
            elif (col['property'] == 'score'):
                if (scanid is not None):
                    col2 = fragal.score
                else:
                    raise ScanRequiredError()
            elif (col['property'] == 'deltappm'):
                if (scanid is not None):
                    col2 = fragal.deltappm
                else:
                    raise ScanRequiredError()
            else:
                col2 = Metabolite.__dict__[col['property']]  #@UndefinedVariable @IgnorePep8
            if (col['direction'] == 'DESC'):
                q = q.order_by(desc(col2))
            elif (col['direction'] == 'ASC'):
                q = q.order_by(asc(col2))

        for r in q[start:(limit + start)]:
            met = r.Metabolite
            row = {
                'metid': met.metid,
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
                'reference': met.reference
            }
            if ('score' in r.keys()):
                row['score'] = r.score
                row['deltappm'] = r.deltappm
            mets.append(row)

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

        for filter in filters:
            if (filter['field'] == 'nhits'
              and filter['comparison'] == 'gt'
              and filter['value'] == 0):
                continue
            if (filter['field'] == 'score'):
                fq = self.extjsgridfilter(fq, Fragment.score, filter)
            elif (filter['field'] == 'deltappm'):
                fq = self.extjsgridfilter(fq, Fragment.deltappm, filter)
            elif (filter['field'] == 'assigned'):
                filter['type'] = 'null'
                fq = fq.join(Peak, and_(Fragment.scanid == Peak.scanid,
                                        Fragment.mz == Peak.mz))
                fq = self.extjsgridfilter(fq, Peak.assigned_metid, filter)
            else:
                fq = fq.join(
                             Metabolite,
                             Fragment.metabolite  #@UndefinedVariable @IgnorePep8
                             )
                fq = self.extjsgridfilter(
                                          fq,
                                          Metabolite.__dict__[filter['field']],  #@UndefinedVariable @IgnorePep8
                                          filter)

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
        """Returns dict with scans key with list of dicts with
         - id
         - rt
         - basepeakintensity
        keys for each lvl1 scan

        and cutoff key with ms_intensity_cutoff
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
            return {
                    'scans': scans,
                    'cutoff': runInfo.ms_intensity_cutoff
                    }
        else:
            return {
                    'scans': scans,
                    'cutoff': None
                    }

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
            cutoff = self.session.query(
                Scan.basepeakintensity * Run.msms_intensity_cutoff
                ).filter(Scan.scanid == scanid).scalar()

        peaks = []
        for peak in self.session.query(Peak).filter_by(scanid=scanid):
            peaks.append({
                'mz': peak.mz,
                'intensity': peak.intensity,
                'assigned_metid': peak.assigned_metid
            })

        return {'peaks': peaks, 'cutoff': cutoff, 'mslevel': scan.mslevel,
                'precursor': {
                    'id': scan.precursorscanid, 'mz': scan.precursormz
                    }
                }

    def fragments(self, scanid, metid, node):
        """Returns fragments of a metabolite on a scan.
        When node is not set then returns metabolites and its lvl2 fragments.
        When node is set then returns children fragments as list
            which have ``node`` as parent fragment.

        Can be used in a Extjs.data.TreeStore if jsonified.

        ``scanid``
            Fragments on scan with this identifier

        ``metid``
            Fragments of metabolite with this identifier

        ``node``
            The fragment identifier to fetch children fragments for.
            Use 'root' for root node.

        Raises FragmentNotFound when no fragment is found
            with the given scanid/metid combination
        """
        def q():
            return self.session.query(Fragment,
                                      Metabolite.mol,
                                      Scan.mslevel).join(Metabolite).join(Scan)

        def fragment2json(row):
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
                'deltappm': frag.deltappm
            }
            if (len(frag.children) > 0):
                f['expanded'] = False
                f['leaf'] = False
            else:
                f['expanded'] = True
                f['leaf'] = True
            return f

        # parent metabolite
        if (node == 'root'):
            structures = []
            for row in q().filter(
                                  Fragment.scanid == scanid
                         ).filter(
                                  Fragment.metid == metid
                         ).filter(
                                  Fragment.parentfragid == 0
                         ):
                structure = fragment2json(row)

                qa = self.session.query(func.count('*')).\
                    filter(Peak.scanid == scanid).\
                    filter(Peak.assigned_metid == metid)
                structure['isAssigned'] = qa.scalar() > 0

                # load children
                structure['children'] = []
                pf = Fragment.parentfragid
                for frow in q().filter(pf == structure['fragid']):
                    structure['expanded'] = True
                    structure['children'].append(fragment2json(frow))
                structures.append(structure)

            if (len(structures) == 0):
                raise FragmentNotFound()
            return {'children': structures, 'expanded': True}
        # fragments
        else:
            fragments = []
            for row in q().filter(Fragment.parentfragid == node):
                fragments.append(fragment2json(row))
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
