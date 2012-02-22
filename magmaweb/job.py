import uuid
import os
import csv
import StringIO
import urllib2
import json

from sqlalchemy import create_engine, and_
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import sessionmaker, aliased
from sqlalchemy.sql import func
from sqlalchemy.sql.expression import desc, asc
from sqlalchemy.orm.exc import NoResultFound
from magmaweb.models import Metabolite, Scan, Peak, Fragment, Run

class ScanRequiredError(Exception):
    """ Raised when a scan identifier is required, but non is supplied"""
    pass

class ScanNotFound(Exception):
    """ Raised when a scan identifier is not found"""
    pass

class FragmentNotFound(Exception):
    """ Raised when a fragment is not found"""
    pass

class JobNotFound(Exception):
    """ Raised when a job with a identifier is not found"""
    def __init__(self, jobid):
        self.jobid = jobid

class JobQuery(object):
    structures = None
    mzxml_filename = None
    mzxml_file = None
    n_reaction_steps = None
    metabolism_types = None
    max_broken_bonds = None
    ionisation_mode = None
    use_fragmentation = None
    use_msms_only = None
    ms_intensity_cutoff = None
    msms_intensity_cutoff = None
    mz_precision = None
    precursor_mz_precision = None
    abs_peak_cutoff = None
    rel_peak_cutoff = None
    max_ms_level = None

    def __eq__(self, other):
        """ Compares all attributes except file objects"""
        return (
                self.structures == other.structures
                and
                self.mzxml_filename == other.mzxml_filename
                and
                self.n_reaction_steps == other.n_reaction_steps
                and
                self.metabolism_types == other.metabolism_types
                and
                self.max_broken_bonds == other.max_broken_bonds
                and
                self.ionisation_mode == other.ionisation_mode
                and
                self.use_fragmentation == other.use_fragmentation
                and
                self.use_msms_only == other.use_msms_only
                and
                self.abs_peak_cutoff == other.abs_peak_cutoff
                and
                self.rel_peak_cutoff == other.rel_peak_cutoff
                and
                self.ms_intensity_cutoff == other.ms_intensity_cutoff
                and
                self.msms_intensity_cutoff == other.msms_intensity_cutoff
                and
                self.mz_precision == other.mz_precision
                and
                self.precursor_mz_precision == other.precursor_mz_precision
                and
                self.max_ms_level == other.max_ms_level
                )

class JobFactory(object):
    """ Factory which can create jobs """
    def __init__(
                 self, jobrootdir, job_script, job_tarball, dbname = 'results.db',
                 submiturl='http://localhost:9998', jobstatefilename = 'job.state'
                 ):
        """
        jobrootdir
            Directory in which jobs are created, retrieved

        dbname
            Sqlite db file name in job directory (default results.db)

        submiturl
            Url of job manager daemon where job can be submitted (default http://localhost:9998)

        job_script
            Local absolute location of script which must be run by job manager daemon

        job_tarball
            Local absolute location of tarball which contains application which job_script will unpack and run

        jobstatefilename
            Filename where job manager daemon writes job state (default job.state)
        """
        self.jobrootdir = jobrootdir
        self.dbname = dbname
        self.submiturl = submiturl
        self.jobstatefilename = jobstatefilename
        self.job_script = job_script
        self.job_tarball = job_tarball

    def fromId(self, jobid):
        """
        Finds job db in job root dir.
        Returns a Job instance

        jobid
            Job identifier

        Raises JobNotFound exception when job is not found by jobid
        """
        engine = create_engine(self.id2url(jobid))
        try:
            engine.connect()
        except OperationalError:
            raise JobNotFound(jobid)
        session = sessionmaker(bind=engine)
        return Job(jobid, session(), self.id2jobdir(jobid))

    def fromDb(self, dbfile):
        """
        A job directory is created and the dbfile is copied into it
        Returns a Job instance

        dbfile
            The sqlite result db

        """
        jobid = uuid.uuid4()

        # create job dir
        os.makedirs(self.id2jobdir(jobid))
        # copy dbfile into job dir
        jobdb = open(self.id2db(jobid), 'wb')
        dbfile.seek(0)
        while 1:
            data = dbfile.read(2 << 16)
            if not data:
                break
            jobdb.write(data)
        jobdb.close()

        return self.fromId(jobid)


    def submitJob2Manager(self, body):
        request = urllib2.Request(
                                  self.submiturl,
                                  json.dumps(body),
                                  { 'Content-Type': 'application/json' }
                                  )
        return urllib2.urlopen(request)

    def submitQuery(self, query):
        """
        A job directory is created, the query input files are written to directory and job is submitted to job manager

        query is JobQuery object

        Returns job identifier
        """

        jobid = uuid.uuid4()
        jobdir = self.id2jobdir(jobid)
        os.makedirs(jobdir)

        # copy mzxml file
        jobmzxml = open(os.path.join(jobdir ,'data.mzxml'), 'wb')
        query.mzxml_file.seek(0)
        while 1:
            data = query.mzxml_file.read(2 << 16)
            if not data:
                break
            jobmzxml.write(data)
        jobmzxml.close()

        # copy metabolites string to file
        metsfile = file(os.path.join(jobdir, 'smiles.txt'), 'w')
        metsfile.write(query.structures)
        metsfile.close()

        body = {
                'jobdir': jobdir+'/',
                'executable': "/bin/sh",
                'prestaged': [
                              self.job_script,
                              self.job_tarball,
                              'data.mzxml', 'smiles.txt'
                              ],
                "poststaged": [self.dbname],
                "stderr": "stderr.txt",
                "stdout": "stdout.txt",
                'arguments': [
                              "magma.sh",
                              "allinone",
                              "--mz_precision", query.mz_precision,
                              "--ms_intensity_cutoff", query.ms_intensity_cutoff,
                              "--msms_intensity_cutoff", query.msms_intensity_cutoff,
                              "--ionisation_mode", query.ionisation_mode,
                              "--n_reaction_steps", query.n_reaction_steps,
                              "--metabolism_types", " ".join(query.metabolism_types),
                              "--max_broken_bonds", query.max_broken_bonds,
                              "--abs_peak_cutoff", query.abs_peak_cutoff,
                              "--rel_peak_cutoff", query.rel_peak_cutoff,
                              "--precursor_mz_precision", query.precursor_mz_precision,
                              'data.mzxml', 'smiles.txt',
                              self.dbname
                              ]
                 }
        if query.use_msms_only:
            body['arguments'].insert(-3,"--use_msms_only")
        if query.use_fragmentation:
            body['arguments'].insert(-3,"--use_fragmentation")

        self.submitJob2Manager(body)

        return jobid

    def state(self, id):
        """ Returns state of job, see ibis org.gridlab.gat.resources.Job JobState enum for possible states."""
        try:
            jobstatefile = open(os.path.join(self.id2jobdir(id), self.jobstatefilename))
            jobstate = jobstatefile.readline().strip()
            jobstatefile.close()
            return jobstate
        except IOError:
            return 'UNKNOWN'

    def id2jobdir(self, id):
        """ Returns job directory based on id and jobrootdir """
        return os.path.join(self.jobrootdir, str(id))

    def id2db(self, id):
        """ Returns sqlite db of job with id """
        return os.path.join(self.id2jobdir(id), self.dbname)

    def id2url(self, id):
        """ Returns sqlalchemy url of sqlite db of job with id """
        # 3rd / is for username:pw@host which sqlite does not need
        return 'sqlite:///' + self.id2db(id)

class Job(object):
    """
    Job contains results database of Magma calculation run
    """
    def __init__(self, id, session, dir):
        """
        jobid
            A UUID of the job
        session
            Sqlalchemy session to read/write to job db
        dir
            Directory where input and output files reside
        """
        self.id = id
        self.session = session
        self.dir = dir

    def maxMSLevel(self):
        """ Returns the maximum nr of MS levels """
        return self.session.query(func.max(Scan.mslevel)).scalar()

    def runInfo(self):
        """ Returns run info"""
        return self.session.query(Run).one()

    def extjsgridfilter(self, q, column, filter):
        """Query helper to convert a extjs grid filter to a sqlalchemy query filter"""
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

    def metabolites(self, start=0, limit=10, sorts=[{"property":"probability", "direction":"DESC"}, {"property":"metid", "direction":"ASC"}], scanid=None, filters=[]):
        """
        Returns dict with total and rows attribute

        start
            Offset
        limit
            Maximum nr of metabolites to return
        scanid
            Only return metabolites that have hits in scan with this identifier. Adds score column.
        filters
            List of dicts which is generated by ExtJS component Ext.ux.grid.FiltersFeature
        sorts
            How to sort metabolites. List of dicts. Eg.
                [{"property":"probability","direction":"DESC"},{"property":"metid","direction":"ASC"}]
        """
        mets = []
        q = self.session.query(Metabolite)

            # custom filters
        fragal = aliased(Fragment)
        if (scanid != None):
            # TODO add score column + order by score
            q = q.add_column(fragal.score).join(fragal.metabolite).filter(
                fragal.parentfragid == 0).filter(fragal.scanid == scanid)

        # add nr_scans column
        stmt = self.session.query(Fragment.metid, func.count('*').label('nr_scans')).filter(
            Fragment.parentfragid == 0).group_by(Fragment.metid).subquery()
        q = q.add_column(stmt.c.nr_scans).outerjoin(stmt, Metabolite.metid == stmt.c.metid)

        for filter in filters:
            # generic filters
            if (filter['field'] == 'nr_scans'):
                col = stmt.c.nr_scans
            elif (filter['field'] == 'score'):
                if (scanid != None):
                    col = fragal.score
                else:
                    raise ScanRequiredError()
            else:
                col = Metabolite.__dict__[filter['field']] #@UndefinedVariable
            q = self.extjsgridfilter(q, col, filter)

        total = q.count()

        for col in sorts:
            if (col['property'] == 'nr_scans'):
                col2 = stmt.c.nr_scans
            elif (col['property'] == 'score'):
                if (scanid != None):
                    col2 = fragal.score
                else:
                    raise ScanRequiredError()
            else:
                col2 = Metabolite.__dict__[col['property']] #@UndefinedVariable
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
                'nr_scans': r.nr_scans,
                'mim': met.mim
            }
            if ('score' in r.keys()):
                row['score'] = r.score
            mets.append(row)

        return { 'total': total, 'rows': mets }

    def metabolites2csv(self, metabolites):
        """
        Converts array of metabolites to csv file handler

        Params:
        metabolites array like metabolites()['rows']

        Return
        StringIO
        """
        csvstr = StringIO.StringIO()
        headers = [
                   'name', 'smiles', 'probability', 'reactionsequence',
                   'nr_scans', 'molformula', 'mim' , 'isquery'
                   ]
        if ('score' in metabolites[0].keys()):
            headers.append('score')

        csvwriter = csv.DictWriter(csvstr, headers, extrasaction='ignore')
        csvwriter.writeheader()
        for m in metabolites:
            m['name'] = m['origin']
            csvwriter.writerow(m)

        return csvstr

    def scansWithMetabolites(self, filters=[], metid=None):
        """Returns id and rt of lvl1 scans which have a fragment in it and for which the filters in params pass

        params:

        metid
            Only return scans that have hits with metabolite with this identifier
        filters
            List which is generated by ExtJS component Ext.ux.grid.FiltersFeature, with columns from Metabolite grid
        """
        fq = self.session.query(Fragment.scanid).filter(Fragment.parentfragid == 0)
        if (metid != None):
            fq = fq.filter(Fragment.metid == metid)
        fq = fq.join(Metabolite)
        for filter in filters:
            if (filter['field'] == 'score'):
                fq = self.extjsgridfilter(fq, Fragment.score, filter)
            elif (filter['field'] != 'nr_scans'):
                fq = self.extjsgridfilter(
                                          fq,
                                          Metabolite.__dict__[filter['field']], #@UndefinedVariable
                                          filter)

        hits = []
        for hit in self.session.query(Scan.rt, Scan.scanid).filter_by(mslevel=1).filter(Scan.scanid.in_(fq)):
            hits.append({
                'id': hit.scanid,
                'rt': hit.rt
            })

        return hits

    def extractedIonChromatogram(self, metid):
        """ Returns extracted ion chromatogram of metabolite with id metid """
        chromatogram = []
        mzq = self.session.query(func.avg(Fragment.mz)).filter(Fragment.metid == metid).filter(Fragment.parentfragid == 0).scalar()
        mzoffset = self.session.query(Run.mz_precision).scalar()
        # fetch max intensity of peaks with mz = mzq+-mzoffset
        for (rt, intens) in self.session.query(Scan.rt, func.max(Peak.intensity)).outerjoin(Peak, and_(Peak.scanid == Scan.scanid, Peak.mz.between(mzq - mzoffset, mzq + mzoffset))).filter(Scan.mslevel == 1).group_by(Scan.rt).order_by(asc(Scan.rt)):
            chromatogram.append({
                'rt': rt,
                'intensity': intens or 0
            })
        return chromatogram

    def chromatogram(self):
        """Returns list of dicts with the id, rt and basepeakintensity for each lvl1 scan"""
        scans = []
        # TODO add left join to find if scan has metabolite hit
        for scan in self.session.query(Scan).filter_by(mslevel=1):
            scans.append({
                'id': scan.scanid,
                'rt': scan.rt,
                'intensity': scan.basepeakintensity
            })

        return scans

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
        if (mslevel != None):
            scanq = scanq.filter(Scan.mslevel == mslevel)

        try:
            scan = scanq.one()
        except NoResultFound:
            raise ScanNotFound()

        # lvl1 scans use absolute cutoff, lvl>1 use ratio of basepeak as cutoff
        if (scan.mslevel == 1):
            cutoff = self.session.query(Run.ms_intensity_cutoff).scalar()
        else:
            cutoff = self.session.query(Scan.basepeakintensity * Run.msms_intensity_cutoff).filter(Scan.scanid == scanid).scalar()

        peaks = []
        for peak in self.session.query(Peak).filter_by(scanid=scanid):
            peaks.append({
                'mz': peak.mz,
                'intensity': peak.intensity
            })

        return { 'peaks': peaks, 'cutoff': cutoff, 'mslevel': scan.mslevel, 'precursor': { 'id': scan.precursorscanid, 'mz': scan.precursormz } }

    def fragments(self, scanid, metid, node=''):
        """Returns dict with metabolites and its lvl2 fragments when node is not set
        When node is set then returns the children fragments as list which have node as parent fragment

        Can be used in a Extjs.data.TreeStore if jsonified.

        scanid
            Fragments on scan with this identifier

        metid
            Fragments of metabolite with this identifier

        node
            The fragment identifier to fetch children fragments for. Optional.

        Raises FragmentNotFound when no fragment is found with scanid/metid combination
        """
        def q():
            return self.session.query(Fragment, Metabolite.mol, Scan.mslevel).join(Metabolite).join(Scan)

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
            }
            if (len(frag.children) > 0):
                f['expanded'] = False
                f['leaf'] = False
            else:
                f['expanded'] = True
                f['leaf'] = True
            return f

        # parent metabolite
        if (node == ''):
            try:
                row = q().filter(
                    Fragment.scanid == scanid).filter(
                    Fragment.metid == metid).filter(
                    Fragment.parentfragid == 0).one()
            except NoResultFound:
                raise FragmentNotFound()
            metabolite = fragment2json(row)
            metabolite['children'] = []
            for frow in q().filter(Fragment.parentfragid == metabolite['fragid']):
                metabolite['expanded'] = True
                metabolite['children'].append(fragment2json(frow))
            return { 'children': metabolite, 'expanded': True}
        # fragments
        else:
            fragments = []
            for row in q().filter(Fragment.parentfragid == node):
                fragments.append(fragment2json(row))
            return fragments

    def stderr(self):
        """Returns stderr text file"""
        return open(os.path.join(self.dir, 'stderr.txt'), 'rb')
