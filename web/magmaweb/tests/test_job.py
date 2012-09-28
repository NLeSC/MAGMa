import unittest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from mock import Mock, patch
from magmaweb.job import JobFactory, Job, make_job_factory, JobQuery
from magmaweb.models import Metabolite, Scan, Peak, Fragment, Run


def initTestingDB(url='sqlite://', dataset='default'):
    """Creates testing db and populates with test data"""
    # default use a in memory db
    engine = create_engine(url)
    session = sessionmaker(bind=engine)
    dbh = session()
    from magmaweb.models import Base
    Base.metadata.create_all(engine)
    if (dataset == 'default'):
        populateTestingDB(dbh)
    elif (dataset == 'useallpeaks'):
        populateWithUseAllPeaks(dbh)
    return dbh


def populateTestingDB(session):
    """Populates test db with data

    Adds 1 metabolite with one fragment.
    Adds 1 metabolite with one fragment which has 2 child fragments
        of which one has another child fragment.
    Run, Scan and Peak are filled to create a working run.

    session
        session connection to db
    """
    session.add(Run(
        n_reaction_steps=2, metabolism_types='phase1,phase2',
        ionisation_mode=-1, skip_fragmentation=True,
        ms_intensity_cutoff=200000.0, msms_intensity_cutoff=0.5,
        mz_precision=10, use_all_peaks=True,
        ms_filename='F123456.mzxml', abs_peak_cutoff=1000,
        rel_peak_cutoff=0.001, max_ms_level=3, precursor_mz_precision=10,
        max_broken_bonds=4, description='My first description'
    ))
    session.add(Metabolite(
        metid=72, mol='Molfile', level=0, probability=1.0,
        reactionsequence='PARENT', smiles='Oc1ccccc1O',
        molformula='C6H6O2', isquery=True, nhits=1,
        origin='pyrocatechol', mim=110.03677, logp=1.231,
        reference='<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=289">CID: 289</a>'
    ))
    session.add(Scan(
        scanid=641, mslevel=1, rt=933.317, lowmz=90.3916, highmz=1197.78,
        basepeakmz=305.034, basepeakintensity=807577.0, totioncurrent=5957620
    ))
    session.add_all([
        # basepeak
        Peak(scanid=641, mz=305.033508300781, intensity=807576.625),
        # peak of metabolite
        Peak(scanid=641, mz=109.0295639038086, intensity=345608.65625)
    ])
    session.add(Fragment(
        fragid=948,
        metid=72,
        scanid=641,
        mz=109.0295639038086,
        mass=110.0367794368,
        score=200,
        parentfragid=0,
        atoms="0,1,2,3,4,5,6,7",
        deltah=-1.0,
        deltappm=-1.84815979523607e-08
    ))
    # fragments of metid=352 + scanid=870
    session.add(Metabolite(
        isquery=True, level=0, metid=352,
        mol="Molfile of dihydroxyphenyl-valerolactone",
        molformula="C11H12O4",
        origin="dihydroxyphenyl-valerolactone",
        probability=1.0, reactionsequence="PARENT",
        smiles="O=C1OC(Cc2ccc(O)c(O)c2)CC1",
        reference='<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=152432">CID: 152432</a>',
        mim=208.07355, logp=2.763, nhits=1
    ))
    session.add_all([Scan(
        scanid=870, mslevel=1, rt=1254.15, lowmz=91.0302, highmz=1171.51,
        basepeakmz=287.023, basepeakintensity=1972180.0, totioncurrent=9265290
    ), Scan(
        scanid=871, mslevel=2, rt=1254.93, lowmz=51.5211, highmz=216.864,
        basepeakmz=163.076, basepeakintensity=279010.0, totioncurrent=809307,
        precursormz=207.0663147, precursorintensity=293096.0,
        precursorscanid=870
    ), Scan(
        scanid=872, mslevel=3, rt=1256.77, lowmz=50.3338, highmz=172.155,
        basepeakmz=119.087, basepeakintensity=17387.0, totioncurrent=236842,
        precursormz=163.0762329, precursorintensity=6163.73,
        precursorscanid=871
    )])
    session.add_all([
        # basepeak
        Peak(scanid=870, mz=287.022979736328, intensity=1972180.625),
        # peak of metabolite
        Peak(scanid=870, mz=207.066284179688, intensity=293095.84375),
        # basepeak and peak of frag 1709
        Peak(scanid=871, mz=163.076232910156, intensity=279010.28125),
        # peak of frag 1708
        Peak(scanid=871, mz=123.04508972168, intensity=211603.046875),
        # basepeak and peak of frag 1710
        Peak(scanid=872, mz=119.086540222168, intensity=17386.958984375),
    ])
    session.add_all([Fragment(
        fragid=1707, metid=352, scanid=870, mass=208.0735588736,
        mz=207.066284179688, score=100, parentfragid=0,
        atoms="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14", deltah=-1,
        deltappm=-8.18675317722029e-09
    ), Fragment(
        fragid=1708, metid=352, scanid=871, mass=123.0446044689,
        mz=123.04508972167969, score=201, parentfragid=1707,
        atoms="6,7,8,9,10,11,12,13,14", deltah=0,
        deltappm=3.943698856902144e-12
    ), Fragment(
        fragid=1709, metid=352, scanid=871, mass=164.08372962939995,
        mz=163.07623291015625, score=65, parentfragid=1707,
        atoms="3,4,5,6,7,8,9,10,11,12,13,14", deltah=-1,
        deltappm=-1.235815738001507e-08
    ), Fragment(
        fragid=1710, scanid=872, metid=352, mass=116.0626002568,
        mz=119.08654022216797, score=4, parentfragid=1709,
        atoms="4,5,6,7,8,9,11,13,14", deltah=3,
        deltappm=5.0781684060061766e-08
    )])

    session.flush()


def populateWithUseAllPeaks(session):
    """ Dataset with multiple fragments of same metabolite on lvl1 scan """
    session.add(Run(
        n_reaction_steps=2, metabolism_types='phase1,phase2',
        ionisation_mode=-1, skip_fragmentation=True,
        ms_intensity_cutoff=200000.0, msms_intensity_cutoff=0.5,
        mz_precision=10, use_all_peaks=False,
        ms_filename='F123456.mzxml', abs_peak_cutoff=1000,
        rel_peak_cutoff=0.001, max_ms_level=3, precursor_mz_precision=10,
        max_broken_bonds=4, description='My second description'
    ))
    session.add(Metabolite(
        metid=12,
        level=1,
        probability=0.119004,
        reactionsequence='sulfation_(aromatic_hydroxyl)',
        smiles='Oc1ccc(CC2OC(=O)CC2)cc1OS(O)(=O)=O',
        molformula='C11H12O7S',
        isquery=False,
        origin='5-(3,4)-dihydroxyphenyl-g-valerolactone (F)',
        mol='Molfile',
        reference='',
        mim=288.0303734299,
        logp=1.9027,
        nhits=1
    ))
    session.add_all([Scan(
        scanid=1,
        mslevel=1,
        rt=0.503165,
        lowmz=286.529,
        highmz=288.239,
        basepeakmz=287.023,
        basepeakintensity=39047000.0,
        totioncurrent=49605000.0,
        precursorscanid=0
    ), Scan(
        scanid=2,
        mslevel=2,
        rt=0.544193333333333,
        lowmz=66.8575,
        highmz=288.026,
        basepeakmz=207.066,
        basepeakintensity=32485600.0,
        totioncurrent=42005700.0,
        precursormz=287.0231323,
        precursorintensity=39047000.0,
        precursorscanid=1
    )])
    session.add_all([Peak(
        scanid=1,
        mz=287.015686035156,
        intensity=1058332.875
    ), Peak(
        scanid=1,
        mz=287.023132324219,
        intensity=39047040.0
    ), Peak(
        scanid=2,
        mz=207.066223144531,
        intensity=32485624.0
    ), Peak(
        scanid=2,
        mz=287.022827148438,
        intensity=6491798.0
    )])
    session.add_all([Fragment(
        fragid=17,
        metid=12,
        scanid=1,
        mz=287.015686035156,
        mass=288.0303734299,
        score=0.0,
        parentfragid=0,
        atoms='0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
        deltah=-1.0,
        deltappm=-7.046696487857745e-09
    ), Fragment(
        fragid=18,
        metid=12,
        scanid=1,
        mz=287.023132324219,
        mass=288.0303734299,
        score=0.5,
        parentfragid=0,
        atoms='0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
        deltah=-1.0,
        deltappm=-7.020570507205176e-09
    ), Fragment(
        fragid=19,
        metid=12,
        scanid=2,
        mz=207.066223144531,
        mass=207.0657338415,
        score=0.5,
        parentfragid=18,
        atoms='0,1,2,3,4,5,6,7,8,9,10,11,12,13,14',
        deltah=0.0,
        deltappm=2.3630267823129625e-12
    ), Fragment(
        fragid=20,
        metid=12,
        scanid=2,
        mz=287.022827148438,
        mass=288.0303734299,
        score=0.0,
        parentfragid=18,
        atoms='0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
        deltah=-1.0,
        deltappm=-7.021641217476223e-09
    )])


class JobFactoryFactoryTestCase(unittest.TestCase):

    def test_norootdir(self):
        with self.assertRaises(TypeError):
            make_job_factory({})

    def test_minimal(self):
        factory = make_job_factory({'jobfactory.root_dir': '/somedir'})

        self.assertEqual(factory.root_dir, '/somedir')
        self.assertEqual(factory.init_script, '')
        self.assertEqual(factory.tarball, None)

    def test_local(self):
        settings = {
                    'jobfactory.root_dir': '/somedir',
                    'jobfactory.init_script': '. /somedir/env/bin/activate'
                    }
        factory = make_job_factory(settings)
        self.assertEqual(factory.root_dir, '/somedir')
        self.assertEqual(factory.init_script, '. /somedir/env/bin/activate')
        self.assertEqual(factory.tarball, None)

    def test_remote(self):
        settings = {
                    'jobfactory.root_dir': '/somedir',
                    'jobfactory.init_script': 'tar -zxf Magma.tar.gz;',
                    'jobfactory.tarball': '/somepath/Magma.tar.gz'
                    }
        factory = make_job_factory(settings)
        self.assertEqual(factory.root_dir, '/somedir')
        self.assertEqual(factory.init_script, 'tar -zxf Magma.tar.gz;')
        self.assertEqual(factory.tarball, '/somepath/Magma.tar.gz')


class JobFactoryTestCase(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.root_dir = tempfile.mkdtemp()
        self.factory = JobFactory(root_dir=self.root_dir)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.root_dir)

    def test_hasrootdir(self):
        self.assertEqual(self.factory.root_dir, self.root_dir)

    def test_id2url(self):
        import uuid
        import os
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jobdbfn = 'sqlite:///'
        jobdbfn += os.path.join(self.root_dir, str(jobid), 'results.db')
        self.assertEqual(self.factory.id2url(jobid), jobdbfn)

    def test_id2jobdir(self):
        import uuid
        import os
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jobdir = os.path.join(self.root_dir, str(jobid))
        self.assertEqual(self.factory.id2jobdir(jobid), jobdir)

    def test_id2db(self):
        import uuid
        import os
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jobdbfn = os.path.join(self.root_dir, str(jobid), 'results.db')
        self.assertEqual(self.factory.id2db(jobid), jobdbfn)

    def test_fromdb(self):
        import uuid
        import os
        dbfile = os.tmpfile()
        dbfile.write('bla')

        job = self.factory.fromDb(dbfile)

        self.assertIsInstance(job, Job)
        self.assertIsInstance(job.id, uuid.UUID)
        self.assertEqual(
            str(job.session.get_bind().url),
            self.factory.id2url(job.id),
            'job has dbsession set to sqlite file in job dir'
        )
        self.assertEqual(
            open(self.factory.id2db(job.id)).read(),
            'bla',
            'query db file content has been copied to job dir'
        )

    def test_fromid(self):
        import uuid
        import os
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        os.mkdir(self.factory.id2jobdir(jobid))
        expected_job_dburl = self.factory.id2url(jobid)
        initTestingDB(expected_job_dburl)

        job = self.factory.fromId(jobid)

        self.assertIsInstance(job, Job)
        self.assertEqual(job.id, jobid)
        self.assertEqual(
            str(job.session.get_bind().url),
            expected_job_dburl,
            'job has dbsession set to sqlite file in job dir'
        )

    def test_fromid_notfound(self):
        import uuid
        jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        from magmaweb.job import JobNotFound
        with self.assertRaises(JobNotFound) as exc:
            self.factory.fromId(jobid)
        self.assertEqual(exc.exception.jobid, jobid)

    def test_submitQuery(self):
        import os
        self.factory.init_script = "# make magma available"
        job = self.factory.fromScratch()

        self.factory.script_fn = 'script.sh'
        cmd = "magma add_structures -t smiles structures.dat results.db\n"
        jobquery = JobQuery(job.id, job.dir, cmd, ['structures.dat'])

        self.factory.submitJob2Manager = Mock()

        jobid = self.factory.submitQuery(jobquery)

        self.assertEqual(jobid, jobquery.id)
        job_script_fn = os.path.join(jobquery.dir, self.factory.script_fn)
        job_script = open(job_script_fn).read()
        exp_script = "# make magma available\n"
        exp_script += "magma add_structures -t smiles "
        exp_script += "structures.dat results.db\n"
        self.assertMultiLineEqual(job_script, exp_script)
        jobmanager_query = {
                            'jobdir': jobquery.dir + '/',
                            'executable': "/bin/sh",
                            'prestaged': [
                                          self.factory.script_fn,
                                          'results.db',
                                          'structures.dat'
                                          ],
                            "poststaged": ['results.db'],
                            "stderr": "stderr.txt",
                            "stdout": "stdout.txt",
                            "time_max": self.factory.time_max,
                            'arguments': [self.factory.script_fn]
                            }
        self.factory.submitJob2Manager.assert_called_with(jobmanager_query)

    def test_submitQuery_with_tarball(self):
        self.factory.tarball = 'Magma-1.1.tar.gz'
        self.factory.submitJob2Manager = Mock()
        job = self.factory.fromScratch()
        jobquery = JobQuery(job.id, job.dir, "", [])

        jobid = self.factory.submitQuery(jobquery)

        self.assertEqual(jobid, jobquery.id)
        jobmanager_query = {
                            'jobdir': jobquery.dir + '/',
                            'executable': "/bin/sh",
                            'prestaged': [
                                          self.factory.script_fn,
                                          'results.db',
                                          # tarball is staged as well
                                          'Magma-1.1.tar.gz'
                                          ],
                            "poststaged": ['results.db'],
                            "stderr": "stderr.txt",
                            "stdout": "stdout.txt",
                            "time_max": self.factory.time_max,
                            'arguments': [self.factory.script_fn]
                            }
        self.factory.submitJob2Manager.assert_called_with(jobmanager_query)

    @patch('urllib2.urlopen')
    def test_submitJob2Manager(self, ua):
        import json
        body = {'foo': 'bar'}

        self.factory.submitJob2Manager(body)

        # TODO replace with ua.assert_called_once_with(<hamcrest matcher>)
        req = ua.call_args[0][0]
        self.assertEquals(req.get_data(), json.dumps(body))
        self.assertEquals(req.get_full_url(), self.factory.submit_url)
        self.assertEquals(req.get_header('Content-type'), 'application/json')

    def test_state(self):
        import uuid
        import os
        jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        jobdir = self.factory.id2jobdir(jobid)
        os.mkdir(jobdir)
        jobstatefile = open(os.path.join(jobdir, 'job.state'), 'w')
        jobstatefile.write('STOPPED')
        jobstatefile.close()

        state = self.factory.state(jobid)

        self.assertEquals(state, 'STOPPED')

    def test_stateNotExisting(self):
        import uuid
        import os
        jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        jobdir = self.factory.id2jobdir(jobid)
        os.mkdir(jobdir)

        state = self.factory.state(jobid)

        self.assertEquals(state, 'UNKNOWN')

    def test_fromScratch(self):
        import uuid

        job = self.factory.fromScratch()

        self.assertIsInstance(job.id, uuid.UUID)
        self.assertEqual(job.maxMSLevel(), 0)

    def test_cloneJob(self):
        job = self.factory.fromScratch()
        newjob = self.factory.cloneJob(job)
        self.assertNotEqual(job.id, newjob.id)


class JobNotFound(unittest.TestCase):
    def test_it(self):
        from magmaweb.job import JobNotFound
        import uuid
        jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        e = JobNotFound(jobid)
        self.assertEqual(e.jobid, jobid)


class JobTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        import tempfile
        import os
        self.jobid = uuid.uuid1()
        # mock job session
        self.session = initTestingDB()
        self.jobdir = tempfile.mkdtemp()
        stderr = open(os.path.join(self.jobdir, 'stderr.txt'), 'w')
        stderr.write('Error log')
        stderr.close()
        self.job = Job(self.jobid, self.session, self.jobdir)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.jobdir)

    def test_construct(self):
        self.assertEqual(self.job.id, self.jobid)
        self.assertEqual(self.job.session, self.session)
        self.assertEqual(self.job.dir, self.jobdir)

    def test_runInfo(self):
        runInfo = self.job.runInfo()
        self.assertEqual(runInfo.n_reaction_steps, 2)
        self.assertEqual(runInfo.metabolism_types, "phase1,phase2")

        self.assertEqual(runInfo.ms_filename, 'F123456.mzxml')
        self.assertEqual(runInfo.abs_peak_cutoff, 1000)
        self.assertEqual(runInfo.rel_peak_cutoff, 0.001)
        self.assertEqual(runInfo.max_ms_level, 3)
        self.assertEqual(runInfo.precursor_mz_precision, 10)

        self.assertEqual(runInfo.ionisation_mode, -1)
        self.assertEqual(runInfo.skip_fragmentation, True)
        self.assertEqual(runInfo.max_broken_bonds, 4)
        self.assertEqual(runInfo.ms_intensity_cutoff, 200000.0)
        self.assertEqual(runInfo.msms_intensity_cutoff, 0.5)
        self.assertEqual(runInfo.mz_precision, 10)
        self.assertEqual(runInfo.use_all_peaks, True)
        self.assertEqual(runInfo.description, 'My first description')

    def test_runInfo_maxrunid(self):
        self.session.add(Run(
            n_reaction_steps=2, metabolism_types='phase1,phase2',
            ionisation_mode=-1, skip_fragmentation=True,
            ms_intensity_cutoff=200000.0, msms_intensity_cutoff=0.5,
            mz_precision=10, use_all_peaks=True,
            ms_filename='F123456.mzxml', abs_peak_cutoff=1000,
            rel_peak_cutoff=0.001, max_ms_level=3, precursor_mz_precision=10,
            max_broken_bonds=4, description='My second description'
        ))

        runInfo = self.job.runInfo()

        self.assertEqual(runInfo.description, 'My second description')

    def test_maxMSLevel(self):
        maxmslevel = self.job.maxMSLevel()
        self.assertEqual(maxmslevel, 3)

    def test_maxMSLevel_withoutscans(self):
        self.session.query(Peak).delete()
        self.session.query(Scan).delete()
        maxmslevel = self.job.maxMSLevel()
        self.assertEqual(maxmslevel, 0)

    def test_extractedionchromatogram(self):
        metid = 72
        eic = self.job.extractedIonChromatogram(metid)
        self.assertEqual(eic, [{
            'rt': 933.317,
            'intensity': 345608.65625
        }, {
            'rt': 1254.15,
            'intensity': 0
        }])

    def test_chromatogram(self):
        expected_chromatogram = {
                                 'scans': [
                                          {'id': 641, 'rt': 933.317,
                                           'intensity': 807577.0, 'ap': 0},
                                          {'id': 870, 'rt': 1254.15,
                                           'intensity': 1972180.0, 'ap': 0}
                                          ],
                                 'cutoff': 200000.0
                                 }
        self.assertEqual(self.job.chromatogram(), expected_chromatogram)

    def test_chromatogram_with_assigned_peaks(self):
        metid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz, metid)

        expected_chromatogram = {
                                 'scans': [
                                          {'id': 641, 'rt': 933.317,
                                           'intensity': 807577.0, 'ap': 1},
                                          {'id': 870, 'rt': 1254.15,
                                           'intensity': 1972180.0, 'ap': 0}
                                          ],
                                 'cutoff': 200000.0
                                 }
        self.assertEqual(self.job.chromatogram(), expected_chromatogram)

    def test_stderr(self):
        log = self.job.stderr()
        self.assertIsInstance(log, file)
        self.assertEqual(log.name, self.jobdir + '/stderr.txt')
        self.assertEqual(log.read(), 'Error log')

    def test_metabolitesTotalCount(self):
        self.assertEqual(self.job.metabolitesTotalCount(), 2)

    def test_jobquery(self):
        jobquery = self.job.jobquery()
        self.assertIsInstance(jobquery, JobQuery)
        self.assertEqual(jobquery.id, self.job.id)
        self.assertEqual(jobquery.dir, self.job.dir)

    def test_set_description(self):
        self.job.description('My second description')
        self.assertEqual(self.job.runInfo().description,
                         'My second description')

    def test_assign_metabolite2peak(self):
        metid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz, metid)

        q = self.session.query(Peak.assigned_metid)
        q = q.filter(Peak.scanid == scanid)
        expected_metid = q.filter(Peak.mz == mz).scalar()
        self.assertEqual(expected_metid, metid)

    def test_assign_metabolite2peak_withoffset(self):
        metid = 72
        scanid = 641
        offset = 8e-7
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz + offset, metid)

        q = self.session.query(Peak.assigned_metid)
        q = q.filter(Peak.scanid == scanid)
        expected_metid = q.filter(Peak.mz == mz).scalar()
        self.assertEqual(expected_metid, metid)

    def test_unassign_metabolite2peak(self):
        metid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz, metid)

        self.job.unassign_metabolite2peak(scanid, mz)

        q = self.session.query(Peak.assigned_metid)
        q = q.filter(Peak.scanid == scanid)
        self.assertIsNone(q.filter(Peak.mz == mz).scalar())

    def test_unassign_metabolite2peak_withoffset(self):
        metid = 72
        scanid = 641
        offset = 8e-7
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz, metid)

        self.job.unassign_metabolite2peak(scanid, mz + offset)

        q = self.session.query(Peak.assigned_metid)
        q = q.filter(Peak.scanid == scanid)
        self.assertIsNone(q.filter(Peak.mz == mz).scalar())


class JobEmptyDatasetTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.jobid = uuid.uuid1()
        # mock job session
        self.session = initTestingDB(dataset=None)
        self.jobdir = '/somedir'
        self.job = Job(self.jobid, self.session, self.jobdir)

    def test_runInfo(self):
        runInfo = self.job.runInfo()
        self.assertIsNone(runInfo)

    def test_maxMSLevel(self):
        maxmslevel = self.job.maxMSLevel()
        self.assertEqual(maxmslevel, 0)

    def test_chromatogram(self):
        self.assertEqual(self.job.chromatogram(), {
                                                   'scans': [],
                                                   'cutoff': None
                                                   })

    def test_stderr(self):
        open_stub = Mock(side_effect=IOError('File not found'))
        with patch('__builtin__.open', open_stub):
            log = self.job.stderr()
            self.assertEqual(log.read(), '')

    def test_metabolitesTotalCount(self):
        self.assertEqual(self.job.metabolitesTotalCount(), 0)

    def test_set_description(self):
        desc = 'My second description'
        self.job.description(desc)
        self.assertEqual(self.job.runInfo().description, desc)


class JobMetabolitesTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB(), '/tmp')

    def test_default(self):
        response = self.job.metabolites()
        self.assertEquals(
            response,
            {
                'total': 2,
                'rows': [{
                    'metid': 72,
                    'isquery': True,
                    'level': 0,
                    'mol': u'Molfile',
                    'molformula': u'C6H6O2',
                    'nhits': None,
                    'nhits': 1,
                    'origin': u'pyrocatechol',
                    'probability': 1.0,
                    'reactionsequence': u'PARENT',
                    'smiles': u'Oc1ccccc1O',
                    'mim': 110.03677, 'logp':1.231,
                    'assigned': False,
                    'reference': '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=289">CID: 289</a>'
                }, {
                    'isquery': True, 'level': 0, 'metid': 352,
                    'mol': u"Molfile of dihydroxyphenyl-valerolactone",
                    'molformula': u"C11H12O4",
                    'nhits': None,
                    'nhits': 1,
                    'origin': u"dihydroxyphenyl-valerolactone",
                    'probability': 1, 'reactionsequence': u"PARENT",
                    'smiles': u"O=C1OC(Cc2ccc(O)c(O)c2)CC1",
                    'mim': 208.07355, 'logp':2.763,
                    'assigned': False,
                    'reference': '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=152432">CID: 152432</a>'
                }]
            }
        )

    def test_scanid(self):
        response = self.job.metabolites(scanid=641)
        self.assertIn('score', response['rows'][0])
        self.assertIn('deltappm', response['rows'][0])
        self.assertEqual(response['total'], 1)

    def test_filteredon_nrscanseq(self):
        response = self.job.metabolites(filters=[{"type": "numeric",
                                                  "comparison": "eq",
                                                  "value": 1,
                                                  "field": "nhits"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscansgt(self):
        response = self.job.metabolites(filters=[{"type": "numeric",
                                                  "comparison": "lt",
                                                  "value": 2,
                                                  "field": "nhits"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscanslt(self):
        response = self.job.metabolites(filters=[{"type": "numeric",
                                                  "comparison": "gt",
                                                  "value": 0,
                                                  "field": "nhits"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_isquery(self):
        response = self.job.metabolites(filters=[{"type": "boolean",
                                                  "value": True,
                                                  "field": "isquery"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_molformula(self):
        response = self.job.metabolites(filters=[{"type": "string",
                                                  "value": "C6",
                                                  "field": "molformula"}])

        self.assertEqual(response['total'], 1)

    def test_filteredon_level(self):
        response = self.job.metabolites(filters=[{"type": "list",
                                                  "value": [0, 1, 2],
                                                  "field": "level"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_score(self):
        filters = [{"type": "numeric", "comparison": "eq",
                    "value": 200, "field": "score"}]

        response = self.job.metabolites(scanid=641, filters=filters)

        self.assertEqual(response['total'], 1)

    def test_filteredon_score_without_scan(self):
        from magmaweb.job import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(filters=[{"type": "numeric",
                                           "comparison": "eq",
                                           "value": 200,
                                           "field": "score"}])

    def test_filteredon_deltappm(self):
        filters = [{"type": "numeric", "comparison": "eq",
                    "value": -1.84815979523607e-08, "field": "deltappm"}]

        response = self.job.metabolites(scanid=641, filters=filters)

        self.assertEqual(response['total'], 1)

    def test_filteredon_deltappm_without_scan(self):
        from magmaweb.job import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(filters=[{"type": "numeric",
                                           "comparison": "eq",
                                           "value": -1.84815979523607e-08,
                                           "field": "deltappm"}])

    def test_filteredon_not_assigned(self):
        response = self.job.metabolites(filters=[{"type": "boolean",
                                                  "value": False,
                                                  "field": "assigned"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_assigned(self):
        response = self.job.metabolites(filters=[{"type": "boolean",
                                                  "value": True,
                                                  "field": "assigned"}])

        self.assertEqual(response['total'], 0)

    def test_sort_probmet(self):
        response = self.job.metabolites(sorts=[{"property": "probability",
                                                "direction": "DESC"},
                                               {"property": "metid",
                                                "direction": "ASC"}])

        self.assertEqual(response['total'], 2)

    def test_sort_nrscans(self):
        response = self.job.metabolites(sorts=[{"property": "nhits",
                                                "direction": "DESC"}])

        self.assertEqual(response['total'], 2)

    def test_sort_assigned(self):
        response = self.job.metabolites(sorts=[{"property": "assigned",
                                                "direction": "DESC"}])

        self.assertEqual(response['total'], 2)

    def test_sort_score(self):
        sorts = [{"property": "score", "direction": "DESC"}]

        response = self.job.metabolites(scanid=641, sorts=sorts)

        self.assertEqual(response['total'], 1)

    def test_sort_score_without_scan(self):
        from magmaweb.job import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(sorts=[{"property": "score",
                                         "direction": "DESC"}])

    def test_sort_deltappm(self):
        sorts = [{"property": "deltappm", "direction": "DESC"}]

        response = self.job.metabolites(scanid=641, sorts=sorts)

        self.assertEqual(response['total'], 1)

    def test_sort_deltappm_without_scan(self):
        sorts = [{"property": "deltappm", "direction": "DESC"}]
        from magmaweb.job import ScanRequiredError

        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(sorts=sorts)


class JobMetabolites2csvTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB(), '/tmp')

    def test_it(self):
        csvfile = self.job.metabolites2csv(self.job.metabolites()['rows'])
        import csv
        import StringIO
        expected_csvfile = StringIO.StringIO()
        csvwriter = csv.DictWriter(expected_csvfile, [
                                                  'name', 'smiles',
                                                  'probability',
                                                  'reactionsequence',
                                                  'nhits', 'molformula',
                                                  'mim',
                                                  'isquery', 'logp',
                                                  'reference'
                                                  ])
        csvwriter.writeheader()
        csvwriter.writerow({
                            'name': 'pyrocatechol', 'smiles': 'Oc1ccccc1O',
                            'probability': 1.0, 'reactionsequence': 'PARENT',
                            'nhits': 1, 'molformula': 'C6H6O2',
                            'isquery': True, 'mim': 110.03677, 'logp': 1.231,
                            'reference': '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=289">CID: 289</a>'
                            })
        csvwriter.writerow({
                            'name': 'dihydroxyphenyl-valerolactone',
                            'smiles': 'O=C1OC(Cc2ccc(O)c(O)c2)CC1',
                            'probability': 1.0, 'reactionsequence': 'PARENT',
                            'nhits': 1, 'molformula': 'C11H12O4',
                            'isquery': True, 'mim': 208.07355, 'logp': 2.763,
                            'reference': '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=152432">CID: 152432</a>',
                            })
        self.assertMultiLineEqual(csvfile.getvalue(),
                                  expected_csvfile.getvalue())

    def test_with_score(self):
        mets = self.job.metabolites(scanid=641)['rows']
        csvfile = self.job.metabolites2csv(mets)
        import csv
        import StringIO
        expected_csvfile = StringIO.StringIO()
        csvwriter = csv.DictWriter(expected_csvfile, [
                                                  'name', 'smiles',
                                                  'probability',
                                                  'reactionsequence',
                                                  'nhits', 'molformula', 'mim',
                                                  'isquery', 'logp',
                                                  'reference', 'score'
                                                  ])
        csvwriter.writeheader()
        csvwriter.writerow({
                            'name': 'pyrocatechol', 'smiles': 'Oc1ccccc1O',
                            'probability': 1.0, 'reactionsequence': 'PARENT',
                            'nhits': 1, 'molformula': 'C6H6O2',
                            'isquery': True, 'score': 200.0, 'mim': 110.03677,
                            'logp': 1.231,
                            'reference': '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=289">CID: 289</a>'
                            })
        self.assertMultiLineEqual(csvfile.getvalue(),
                                  expected_csvfile.getvalue())


class JobMetabolites2sdfTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB(), '/tmp')

    def test_it(self):
        sdffile = self.job.metabolites2sdf(self.job.metabolites()['rows'])

        expected_sdf = """Molfile> <name>
pyrocatechol

> <smiles>
Oc1ccccc1O

> <probability>
1.0

> <reactionsequence>
PARENT

> <nhits>
1

> <molformula>
C6H6O2

> <mim>
110.03677

> <logp>
1.231

> <reference>
<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=289">CID: 289</a>

$$$$
Molfile of dihydroxyphenyl-valerolactone> <name>
dihydroxyphenyl-valerolactone

> <smiles>
O=C1OC(Cc2ccc(O)c(O)c2)CC1

> <probability>
1.0

> <reactionsequence>
PARENT

> <nhits>
1

> <molformula>
C11H12O4

> <mim>
208.07355

> <logp>
2.763

> <reference>
<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=152432">CID: 152432</a>

$$$$
"""

        self.assertMultiLineEqual(sdffile, expected_sdf)

    def test_with_scan(self):
        """Include score prop"""
        mets = self.job.metabolites(scanid=641)['rows']
        sdffile = self.job.metabolites2sdf(mets)

        expected_sdf = """Molfile> <name>
pyrocatechol

> <smiles>
Oc1ccccc1O

> <probability>
1.0

> <reactionsequence>
PARENT

> <nhits>
1

> <molformula>
C6H6O2

> <mim>
110.03677

> <logp>
1.231

> <reference>
<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=289">CID: 289</a>

> <score>
200.0

$$$$
"""

        self.assertMultiLineEqual(sdffile, expected_sdf)


class JobScansWithMetabolitesTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB(), '/tmp')

    def test_metid(self):
        response = self.job.scansWithMetabolites(metid=72)
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_all(self):
        response = self.job.scansWithMetabolites()
        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_molformula(self):
        filters = [{"type": "string", "value": "C6", "field": "molformula"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_nrscans(self):
        filters = [{"type": "numeric", "value": "1",
                    "comparison": "eq", "field": "nhits"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_score(self):
        filters = [{"type": "numeric", "value": "200",
                    "comparison": "eq", "field": "score"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317}
        ])

    def test_filteredon_not_assigned(self):
        filters = [{"type": "boolean", "value": False, "field": "assigned"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_filteredon_assigned(self):
        filters = [{"type": "boolean", "value": True, "field": "assigned"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [])

    def test_filteredon_nrscansgt_and_molformula(self):
        filters = [{"type": "numeric", "comparison": "lt",
                    "value": 2, "field": "nhits"},
                   {"type": "string", "value": "C6", "field": "molformula"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(len(response), 1)


class JobMSpectraTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB(), '/tmp')

    def test_scanonly(self):
        self.assertEqual(
            self.job.mspectra(641),
            {
                'peaks': [
                    {'intensity': 345608.65625,
                     'mz': 109.0295639038086,
                     'assigned_metid': None},
                    {'intensity': 807576.625,
                     'mz': 305.033508300781,
                     'assigned_metid': None}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': {'id': None, 'mz': None}
            }
        )

    def test_withmslevel(self):
        self.assertEqual(
            self.job.mspectra(641, 1),
            {
                'peaks': [
                    {'intensity': 345608.65625,
                     'mz': 109.0295639038086,
                     'assigned_metid': None},
                    {'intensity': 807576.625,
                     'mz': 305.033508300781,
                     'assigned_metid': None}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': {'id': None, 'mz': None}
            }
        )

    def test_notfound(self):
        from magmaweb.job import ScanNotFound
        with self.assertRaises(ScanNotFound):
            self.job.mspectra(123)

    def test_lvl2scan(self):
        response = self.job.mspectra(871, 2)
        self.assertEqual(response, {
            'peaks': [
                {'intensity': 211603.046875,
                 'mz': 123.04508972168,
                 'assigned_metid': None},
                {'intensity': 279010.28125,
                 'mz': 163.076232910156,
                 'assigned_metid': None}
            ],
            'cutoff': 139505.0,
            'mslevel': 2,
            'precursor': {'id': 870, 'mz': 207.0663147}
        })

    def test_withassigned_met2peak(self):
        self.job.assign_metabolite2peak(641, 109.0295639038086, 72)
        self.assertEqual(
            self.job.mspectra(641),
            {
                'peaks': [
                    {'intensity': 345608.65625,
                     'mz': 109.0295639038086,
                     'assigned_metid': 72},
                    {'intensity': 807576.625,
                     'mz': 305.033508300781,
                     'assigned_metid': None}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': {'id': None, 'mz': None}
            }
        )


class JobFragmentsTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB(), '/tmp')

    def test_metabolitewithoutfragments(self):
        self.maxDiff = None
        response = self.job.fragments(metid=72, scanid=641, node='root')
        self.assertEqual(response, {
            'children': [{
                'atoms': u'0,1,2,3,4,5,6,7',
                'children': [],
                'deltah': -1.0,
                'deltappm': -1.84815979523607e-08,
                'expanded': True,
                'fragid': 948,
                'leaf': True,
                'mass': 110.0367794368,
                'metid': 72,
                'mol': u'Molfile',
                'mslevel': 1,
                'mz': 109.0295639038086,
                'scanid': 641,
                'score': 200.0,
                'isAssigned': False,
            }], 'expanded': True
        })

    def test_metabolitewithfragments(self):
        self.maxDiff = None
        response = self.job.fragments(metid=352, scanid=870, node='root')
        self.assertEqual(response, {
            'children': [{
                'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14',
                'children': [{
                    'atoms': "6,7,8,9,10,11,12,13,14",
                    'deltah': 0,
                    'deltappm': 3.943698856902144e-12,
                    'expanded': True,
                    'fragid': 1708,
                    'leaf': True,
                    'mass': 123.0446044689,
                    'metid': 352,
                    'mol': "Molfile of dihydroxyphenyl-valerolactone",
                    'mslevel': 2,
                    'mz': 123.04508972167969,
                    'scanid': 871,
                    'score': 201,
                }, {
                    'atoms': "3,4,5,6,7,8,9,10,11,12,13,14",
                    'deltah': -1,
                    'deltappm': -1.235815738001507e-08,
                    'expanded': False,
                    'fragid': 1709,
                    'leaf': False,
                    'mass': 164.08372962939995,
                    'metid': 352,
                    'mol': "Molfile of dihydroxyphenyl-valerolactone",
                    'mslevel': 2,
                    'mz': 163.07623291015625,
                    'scanid': 871,
                    'score': 65
                }],
                'deltah': -1,
                'deltappm': -8.18675317722029e-09,
                'expanded': True,
                'fragid': 1707,
                'leaf': False,
                'mass': 208.0735588736,
                'metid': 352,
                'mol': u'Molfile of dihydroxyphenyl-valerolactone',
                'mslevel': 1,
                'mz': 207.066284179688,
                'scanid': 870,
                'score': 100,
                'isAssigned': False
            }], 'expanded': True
        })

    def test_lvl3fragments(self):
        response = self.job.fragments(metid=352, scanid=870, node=1709)
        self.assertEqual(response, [{
            'atoms': "4,5,6,7,8,9,11,13,14",
            'deltah': 3,
            'expanded': True,
            'fragid': 1710,
            'leaf': True,
            'mass': 116.0626002568,
            'metid': 352,
            'mol': "Molfile of dihydroxyphenyl-valerolactone",
            'mslevel': 3,
            'mz': 119.08654022216797,
            'scanid': 872,
            'score': 4,
            'deltappm': 5.0781684060061766e-08
        }])

    def test_metabolitewithassignedpeak(self):
        self.job.assign_metabolite2peak(641, 109.0295639038086, 72)
        response = self.job.fragments(metid=72, scanid=641, node='root')
        self.assertEqual(response, {
            'children': [{
                'atoms': u'0,1,2,3,4,5,6,7',
                'children': [],
                'deltah': -1.0,
                'deltappm': -1.84815979523607e-08,
                'expanded': True,
                'fragid': 948,
                'leaf': True,
                'mass': 110.0367794368,
                'metid': 72,
                'mol': u'Molfile',
                'mslevel': 1,
                'mz': 109.0295639038086,
                'scanid': 641,
                'score': 200.0,
                'isAssigned': True,
            }], 'expanded': True
        })

    def test_badfragment(self):
        from magmaweb.job import FragmentNotFound
        with self.assertRaises(FragmentNotFound):
            self.job.fragments(metid=70002, scanid=641, node='root')


class JobWithAllPeaksTestCase(unittest.TestCase):

    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(),
                       initTestingDB(dataset='useallpeaks'),
                       '/tmp')

    def test_default(self):
        response = self.job.metabolites()
        self.assertEquals(
            response,
            {
                'total': 1,
                'rows': [{
                    'metid': 12,
                    'isquery': False,
                    'level': 1,
                    'mol': u'Molfile',
                    'molformula': u'C11H12O7S',
                    'nhits': None,
                    'nhits': 1,
                    'origin': u'5-(3,4)-dihydroxyphenyl-g-valerolactone (F)',
                    'probability': 0.119004,
                    'reactionsequence': u'sulfation_(aromatic_hydroxyl)',
                    'smiles': u'Oc1ccc(CC2OC(=O)CC2)cc1OS(O)(=O)=O',
                    'mim': 288.0303734299, 'logp':1.9027,
                    'assigned': False,
                    'reference': ''
                }]
            }
        )

    def test_lvl1fragments(self):
        response = self.job.fragments(metid=12, scanid=1, node='root')
        self.assertEqual(response, {
            'children': [{
                'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
                'children': [],
                'deltah': -1.0,
                'deltappm': -7.046696487857745e-09,
                'expanded': True,
                'fragid': 17,
                'leaf': True,
                'mass': 288.0303734299,
                'metid': 12,
                'mol': u'Molfile',
                'mslevel': 1,
                'mz': 287.015686035156,
                'scanid': 1,
                'score': 0.0,
                'isAssigned': False
            }, {
                'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
                'deltah': -1.0,
                'deltappm': -7.020570507205176e-09,
                'expanded': True,
                'fragid': 18,
                'leaf': False,
                'mz': 287.023132324219,
                'metid': 12,
                'mol': u'Molfile',
                'mslevel': 1,
                'mass': 288.0303734299,
                'scanid': 1,
                'score': 0.5,
                'isAssigned': False,
                'children': [{
                    'fragid': 19,
                    'metid': 12,
                    'scanid': 2,
                    'mz': 207.066223144531,
                    'mass': 207.0657338415,
                    'score': 0.5,
                    'mol': u'Molfile',
                    'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14',
                    'expanded': True,
                    'leaf': True,
                    'mslevel': 2,
                    'deltah': 0.0,
                    'deltappm': 2.3630267823129625e-12
                }, {
                    'fragid': 20,
                    'metid': 12,
                    'scanid': 2,
                    'mz': 287.022827148438,
                    'mass': 288.0303734299,
                    'score': 0.0,
                    'mol': u'Molfile',
                    'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
                    'expanded': True,
                    'leaf': True,
                    'mslevel': 2,
                    'deltah': -1.0,
                    'deltappm': -7.021641217476223e-09
                }]
            }], 'expanded': True
        })


class JobQueryTestCase(unittest.TestCase):
    def setUp(self):
        self.jobid = 12345
        self.jobdir = '/somedir'
        self.jobquery = JobQuery(id=self.jobid, dir=self.jobdir)

    def test_eq(self):
        self.assertEqual(JobQuery(id=self.jobid, dir=self.jobdir),
                         self.jobquery)

    def test_eq_id(self):
        self.assertNotEqual(JobQuery(id=67890, dir=self.jobdir),
                            self.jobquery)

    def test_eq_dir(self):
        self.assertNotEqual(JobQuery(id=self.jobid, dir='/otherdir'),
                            self.jobquery)

    def test_eq_script(self):
        job1 = JobQuery(id=self.jobid, dir=self.jobdir, script='b')
        self.assertNotEqual(job1, self.jobquery)

    def test_eq_prestaged(self):
        job1 = JobQuery(id=self.jobid, dir=self.jobdir, prestaged=[1])
        self.assertNotEqual(job1, self.jobquery)

    def test_repr(self):
        jq = JobQuery('x', 'y', 'z', [123])
        self.assertEqual(jq.__repr__(), "JobQuery('x', 'y', 'z', [123])")

    def test_escape_single_quote(self):
        jq = JobQuery('x', '/y')
        self.assertEquals(jq.escape("'"), '&#39;')


class JobQueryFileTestCase(unittest.TestCase):
    def test_valid(self):
        from cgi import FieldStorage
        f = FieldStorage()
        df = JobQuery.File().deserialize(None, f)

        self.assertEquals(f, df)

    def test_null(self):
        from colander import null
        self.assertEquals(JobQuery.File().deserialize(None, null), null)

    def test_emptystring(self):
        from colander import null
        self.assertEquals(JobQuery.File().deserialize(None, ''), null)

    def test_invalid(self):
        from colander import Invalid, SchemaNode
        n = SchemaNode(JobQuery.File(), name='filefield')
        with self.assertRaises(Invalid) as e:
            n.deserialize(12345)

        self.assertDictEqual(e.exception.asdict(),
                             {'filefield': '12345 is not a cgi.FieldStorage'})

    def test_serialize(self):
        self.assertEquals(JobQuery.File().serialize(None, 12345), 12345)


class JobQueryActionTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        import tempfile
        self.jobid = uuid.uuid1()
        self.jobdir = tempfile.mkdtemp()
        self.jobquery = JobQuery(id=self.jobid, dir=self.jobdir, script='')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.jobdir)

    def fetch_file(self, filename):
        import os.path
        return file(os.path.join(self.jobdir, filename)).read()


class JobQueryAddStructuresTestCase(JobQueryActionTestCase):
    def test_structures_as_string(self):
        params = {'structure_format': 'smiles', 'structures': 'CCO Ethanol'}
        query = self.jobquery.add_structures(params)

        sf = 'structures.dat'
        script = "{magma} add_structures -t 'smiles' structures.dat {db}\n"
        expected_query = JobQuery(
                          id=self.jobid,
                          dir=self.jobdir,
                          prestaged=[sf],
                          script=script
                          )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(params['structures'], self.fetch_file(sf))

    def test_structures_as_file(self):
        import tempfile
        from cgi import FieldStorage
        sfile = tempfile.TemporaryFile()
        sfile.write('foo')
        sfile.flush()
        sfield = FieldStorage()
        sfield.file = sfile
        params = {'structure_format': 'smiles', 'structures_file': sfield}

        query = self.jobquery.add_structures(params)

        script = "{magma} add_structures -t 'smiles' structures.dat {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': ['structures.dat'],
                          'script': script
                          })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('structures.dat'))

    def test_with_metabolize(self):
        from webob.multidict import MultiDict
        params = MultiDict(
                           structure_format='smiles',
                           structures='CCO Ethanol',
                           metabolize='on',
                           n_reaction_steps=2,
                           metabolism_types='phase1'
                           )
        params.add('metabolism_types', 'phase2')

        query = self.jobquery.add_structures(params)

        sf = 'structures.dat'
        script = "{magma} add_structures -t 'smiles' structures.dat {db} |"
        script += "{magma} metabolize -s '2' -m 'phase1,phase2' -j - {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [sf],
                          'script': script
                          })
        self.assertEqual(query, expected_query)

    def test_with_annotate(self):
        params = {
                  'structure_format': 'smiles',
                  'structures': 'CCO Ethanol',
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 0.001,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 0.1,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4
                  }
        query = self.jobquery.add_structures(params, True)

        sf = 'structures.dat'
        script = "{magma} add_structures -t 'smiles' structures.dat {db} |"
        script += "{magma} annotate -p '0.001' -c '200000.0' -d '0.1' -i '1'"
        script += " -b '4' --precursor_mz_precision '0.005' -j - {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [sf],
                          'script': script
                          })
        self.assertEqual(query, expected_query)

    def test_with_metabolize_and_annotate(self):
        from webob.multidict import MultiDict
        params = MultiDict(
                           structure_format='smiles',
                           structures='CCO Ethanol',
                           metabolize='on',
                           n_reaction_steps=2,
                           metabolism_types='phase2',
                           precursor_mz_precision=0.005,
                           mz_precision=0.001,
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=0.1,
                           ionisation_mode=1,
                           max_broken_bonds=4
                           )
        query = self.jobquery.add_structures(params, True)

        sf = 'structures.dat'
        script = "{magma} add_structures -t 'smiles' structures.dat {db} |"
        script += "{magma} metabolize -s '2' -m 'phase2' -j - {db} |"
        script += "{magma} annotate -p '0.001' -c '200000.0' -d '0.1' -i '1'"
        script += " -b '4' --precursor_mz_precision '0.005' -j - {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [sf],
                          'script': script
                          })
        self.assertEqual(query, expected_query)

    def test_without_structures(self):
        params = {
                  'structure_format': 'smiles',
                  'structures_file': '',
                  'structures': ''
                  }
        from colander import Invalid
        with self.assertRaises(Invalid) as e:
            self.jobquery.add_structures(params, False)

        s = 'Either structures or structure_file must be set'
        sf = 'Either structures or structure_file must be set'
        expected = {'structures': s, 'structures_file': sf}
        self.assertDictEqual(e.exception.asdict(), expected)

    def test_with_string_and_file(self):
        import tempfile
        from cgi import FieldStorage
        sfile = tempfile.TemporaryFile()
        sfile.write('foo')
        sfile.flush()
        sfield = FieldStorage()
        sfield.file = sfile
        params = {
                  'structure_format': 'smiles',
                  'structures_file': sfield,
                  'structures': 'bar'
        }

        self.jobquery.add_structures(params)
        # File is kept and string is ignored
        self.assertMultiLineEqual('foo', self.fetch_file('structures.dat'))


class JobQueryAddMSDataTestCase(JobQueryActionTestCase):

    def test_it(self):
        import tempfile
        from cgi import FieldStorage
        msfile = tempfile.TemporaryFile()
        msfile.write('foo')
        msfile.flush()
        msfield = FieldStorage()
        msfield.file = msfile
        params = {
                  'ms_data_format': 'mzxml',
                  'ms_data_file': msfield,
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'rel_peak_cutoff': 0.01
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -l '3' -a '1000.0' -r '0.01' ms_data.dat {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': ['ms_data.dat'],
                          'script': script
                          })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_with_annotate(self):
        import tempfile
        from cgi import FieldStorage
        msfile = tempfile.TemporaryFile()
        msfile.write('foo')
        msfile.flush()
        msfield = FieldStorage()
        msfield.file = msfile
        params = {
                  'ms_data_format': 'mzxml',
                  'ms_data_file': msfield,
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'rel_peak_cutoff': 0.01,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 0.001,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 0.1,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4
                  }

        query = self.jobquery.add_ms_data(params, True)

        script = "{magma} read_ms_data --ms_data_format 'mzxml' "
        script += "-l '3' -a '1000.0' -r '0.01' ms_data.dat {db}\n"
        script += "{magma} annotate -p '0.001' -c '200000.0' -d '0.1'"
        script += " -i '1' -b '4' --precursor_mz_precision '0.005' {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': ['ms_data.dat'],
                          'script': script
                          })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))


class JobQueryMetabolizeTestCase(JobQueryActionTestCase):

    def test_it(self):
        from webob.multidict import MultiDict
        params = MultiDict(
                           n_reaction_steps=2,
                           metabolism_types='phase1'
                           )

        query = self.jobquery.metabolize(params)

        script = "{magma} metabolize -s '2' -m 'phase1' {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [],
                          'script': script
                          })
        self.assertEqual(query, expected_query)

    def test_with_annotate(self):
        from webob.multidict import MultiDict
        params = MultiDict([
                            ('n_reaction_steps', 2),
                            ('metabolism_types', 'phase1'),
                            ('metabolism_types', 'phase2'),
                            ('precursor_mz_precision', 0.005),
                            ('mz_precision', 0.001),
                            ('ms_intensity_cutoff', 200000),
                            ('msms_intensity_cutoff', 0.1),
                            ('ionisation_mode', 1),
                            ('max_broken_bonds', 4)
                            ])
        query = self.jobquery.metabolize(params, True)

        script = "{magma} metabolize -s '2' -m 'phase1,phase2' {db} |"
        script += "{magma} annotate -p '0.001' -c '200000.0' -d '0.1'"
        script += " -i '1' -b '4' --precursor_mz_precision '0.005' -j - {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [],
                          'script': script
                          })
        self.assertEqual(query, expected_query)


class JobQueryMetabolizeOneTestCase(JobQueryActionTestCase):

    def test_it(self):
        from webob.multidict import MultiDict
        params = MultiDict(
                  metid=123,
                  n_reaction_steps=2,
                  metabolism_types='phase1'
                  )
        params.add('metabolism_types', 'phase2')

        query = self.jobquery.metabolize_one(params)

        script = "echo '123' | {magma} metabolize -j - -s '2'"
        script += " -m 'phase1,phase2' {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [],
                          'script': script
                          })
        self.assertEqual(query, expected_query)

    def test_with_annotate(self):
        from webob.multidict import MultiDict
        params = MultiDict(
                  metid=123,
                  n_reaction_steps=2,
                  metabolism_types='phase1',
                  precursor_mz_precision=0.005,
                  mz_precision=0.001,
                  ms_intensity_cutoff=200000,
                  msms_intensity_cutoff=0.1,
                  ionisation_mode=1,
                  max_broken_bonds=4
                  )

        query = self.jobquery.metabolize_one(params, True)

        script = "echo '123' | {magma} metabolize -j - -s '2' "
        script += "-m 'phase1' {db} |"
        script += "{magma} annotate -p '0.001' -c '200000.0' -d '0.1'"
        script += " -i '1' -b '4' --precursor_mz_precision '0.005' -j - {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [],
                          'script': script
                          })
        self.assertEqual(query, expected_query)


class JobQueryAnnotateTestCase(JobQueryActionTestCase):

    def test_it(self):
        params = {
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 0.001,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 0.1,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4
                  }

        query = self.jobquery.annotate(params)

        script = "{magma} annotate -p '0.001' -c '200000.0' -d '0.1' -i '1'"
        script += " -b '4' --precursor_mz_precision '0.005' {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [],
                          'script': script
                          })
        self.assertEqual(query, expected_query)

    def test_all_peaks_skip(self):
        params = {
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 0.001,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 0.1,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4,
                  'use_all_peaks': 'on'
                  }

        query = self.jobquery.annotate(params)

        script = "{magma} annotate -p '0.001' -c '200000.0' -d '0.1'"
        script += " -i '1' -b '4' --precursor_mz_precision '0.005' -u {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [],
                          'script': script
                          })
        self.assertEqual(query, expected_query)

    def test_no_fragmentation(self):
        params = {
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 0.001,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 0.1,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4,
                  'skip_fragmentation': 'on'
                  }

        query = self.jobquery.annotate(params)

        script = "{magma} annotate -p '0.001' -c '200000.0' -d '0.1'"
        script += " -i '1' -b '4' --precursor_mz_precision '0.005' -f {db}\n"
        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': [],
                          'script': script
                          })
        self.assertEqual(query, expected_query)


class JobQueryAllInOneTestCase(JobQueryActionTestCase):

    def test_it(self):
        self.maxDiff = 100000
        import tempfile
        from cgi import FieldStorage
        ms_data_file = tempfile.NamedTemporaryFile()
        ms_data_file.write('foo')
        ms_data_file.flush()
        msfield = FieldStorage()
        msfield.file = ms_data_file
        from webob.multidict import MultiDict
        params = MultiDict(
                           n_reaction_steps=2,
                           ionisation_mode=1,
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=0.1,
                           abs_peak_cutoff=1000,
                           rel_peak_cutoff=0.01,
                           precursor_mz_precision=0.005,
                           max_broken_bonds=4,
                           mz_precision=0.001,
                           metabolism_types='phase1',
                           max_ms_level=3,
                           structures='C1CCCC1 comp1',
                           ms_data_file=msfield,
                           structure_format='smiles',
                           ms_data_format='mzxml'
                           )
        params.add('metabolism_types', 'phase2')

        query = self.jobquery.allinone(params)

        expected_script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        expected_script += " -l '3' -a '1000.0' -r '0.01' ms_data.dat {db}\n"

        expected_script += "{magma} add_structures -t 'smiles'"
        expected_script += " structures.dat {db}\n"

        expected_script += "{magma} metabolize -s '2' -m 'phase1,phase2'"
        expected_script += " {db}\n"

        expected_script += "{magma} annotate -p '0.001' -c '200000.0' -d '0.1'"
        expected_script += " -i '1' -b '4' --precursor_mz_precision '0.005'"
        expected_script += " {db}\n"

        expected_query = JobQuery(**{
                          'id': self.jobid,
                          'dir': self.jobdir,
                          'prestaged': ['ms_data.dat', 'structures.dat'],
                          'script': expected_script
                          })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(params['structures'],
                                  self.fetch_file('structures.dat'))
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))
