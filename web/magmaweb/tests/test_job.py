import datetime
import os
import uuid
import unittest
from mock import Mock, patch
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import transaction
from magmaweb.job import JobFactory, Job, JobDb, make_job_factory, JobQuery
from magmaweb.models import Metabolite, Scan, Peak, Fragment, Run
import magmaweb.user as mu


def initTestingDB(url='sqlite://', dataset='default'):
    """Creates testing db and populates with test data"""
    # default use a in memory db
    engine = create_engine(url)
    session = sessionmaker(bind=engine)
    dbh = session()
    from magmaweb.models import Base
    Base.metadata.create_all(engine)  # @UndefinedVariable
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
    url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
    url1 += '?cid=289">CID: 289</a>'
    url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
    url2 += '?cid=152432">CID: 152432</a>'
    session.add(Run(
        n_reaction_steps=2, metabolism_types='phase1,phase2',
        ionisation_mode=-1, skip_fragmentation=True,
        ms_intensity_cutoff=200000.0, msms_intensity_cutoff=50,
        mz_precision=10, mz_precision_abs=0.002, use_all_peaks=True,
        ms_filename='F123456.mzxml', abs_peak_cutoff=1000,
        max_ms_level=3, precursor_mz_precision=10,
        max_broken_bonds=4, description='My first description',
        max_water_losses=1,
    ))
    session.add(Metabolite(
        metid=72, mol='Molfile', level=0, probability=1.0,
        reactionsequence=['PARENT'], smiles='Oc1ccccc1O',
        molformula='C6H6O2', isquery=True, nhits=1,
        origin='pyrocatechol', mim=110.03677, logp=1.231,
        reference=url1
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
        deltappm=-1.84815979523607e-08,
        formula="C5H4",
    ))
    # fragments of metid=352 + scanid=870
    session.add(Metabolite(
        isquery=True, level=0, metid=352,
        mol="Molfile of dihydroxyphenyl-valerolactone",
        molformula="C11H12O4",
        origin="dihydroxyphenyl-valerolactone",
        probability=1.0, reactionsequence="PARENT\nCHILD\n",
        smiles="O=C1OC(Cc2ccc(O)c(O)c2)CC1",
        reference=url2,
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
        deltappm=-8.18675317722029e-09,
        formula="C3H5O3",
    ), Fragment(
        fragid=1708, metid=352, scanid=871, mass=123.0446044689,
        mz=123.04508972167969, score=201, parentfragid=1707,
        atoms="6,7,8,9,10,11,12,13,14", deltah=0,
        deltappm=3.943698856902144e-12,
        formula="C3H5O3",
    ), Fragment(
        fragid=1709, metid=352, scanid=871, mass=164.08372962939995,
        mz=163.07623291015625, score=65, parentfragid=1707,
        atoms="3,4,5,6,7,8,9,10,11,12,13,14", deltah=-1,
        deltappm=-1.235815738001507e-08,
        formula="C3H5O3",
    ), Fragment(
        fragid=1710, scanid=872, metid=352, mass=116.0626002568,
        mz=119.08654022216797, score=4, parentfragid=1709,
        atoms="4,5,6,7,8,9,11,13,14", deltah=3,
        deltappm=5.0781684060061766e-08,
        formula="C3H5O3",
    )])

    session.flush()


def populateWithUseAllPeaks(session):
    """ Dataset with multiple fragments of same metabolite on lvl1 scan """
    session.add(Run(
        n_reaction_steps=2, metabolism_types='phase1,phase2',
        ionisation_mode=-1, skip_fragmentation=True,
        ms_intensity_cutoff=200000.0, msms_intensity_cutoff=50,
        mz_precision=10, mz_precision_abs=0.002, use_all_peaks=False,
        ms_filename='F123456.mzxml', abs_peak_cutoff=1000,
        max_ms_level=3, precursor_mz_precision=10,
        max_broken_bonds=4, description='My second description',
        max_water_losses=1,
    ))
    session.add(Metabolite(
        metid=12,
        level=1,
        probability=0.119004,
        reactionsequence=['sulfation_(aromatic_hydroxyl)'],
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
        deltappm=-7.046696487857745e-09,
        formula="C4H6O2",
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
        deltappm=-7.020570507205176e-09,
        formula="C4H6O2",
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
        deltappm=2.3630267823129625e-12,
        formula="C4H6O2",
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
        deltappm=-7.021641217476223e-09,
        formula="C4H6O2",
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
        settings = {'jobfactory.root_dir': '/somedir',
                    'jobfactory.init_script': '. /somedir/env/bin/activate'
                    }
        factory = make_job_factory(settings)
        self.assertEqual(factory.root_dir, '/somedir')
        self.assertEqual(factory.init_script, '. /somedir/env/bin/activate')
        self.assertEqual(factory.tarball, None)

    def test_remote(self):
        settings = {'jobfactory.root_dir': '/somedir',
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

        # fill user db
        transaction.begin()
        engine = create_engine('sqlite:///:memory:')
        mu.DBSession.configure(bind=engine)
        mu.Base.metadata.create_all(engine)  # @UndefinedVariable
        mu.DBSession().add(mu.User('bob', 'Bob', 'bob@example.com'))
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        mu.DBSession().add(mu.JobMeta(jobid, owner='bob'))

    def tearDown(self):
        import shutil
        shutil.rmtree(self.root_dir)
        mu.DBSession.remove()

    def test_hasrootdir(self):
        self.assertEqual(self.factory.root_dir, self.root_dir)

    def test_id2url(self):
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jobdbfn = 'sqlite:///'
        jobdbfn += os.path.join(self.root_dir, str(jobid), 'results.db')
        self.assertEqual(self.factory.id2url(jobid), jobdbfn)

    def test_id2jobdir(self):
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jobdir = os.path.join(self.root_dir, str(jobid))
        self.assertEqual(self.factory.id2jobdir(jobid), jobdir)

    def test_id2db(self):
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jobdbfn = os.path.join(self.root_dir, str(jobid), 'results.db')
        self.assertEqual(self.factory.id2db(jobid), jobdbfn)

    def test_fromdb(self):
        # mock/stub private methods which do external calls
        self.factory._makeJobDir = Mock(return_value='/mydir')
        self.factory._copyFile = Mock()
        self.factory._makeJobSession = Mock(return_value=initTestingDB())

        dbfile = os.tmpfile()

        job = self.factory.fromDb(dbfile, 'bob')

        self.assertIsInstance(job.id, uuid.UUID)
        self.assertEqual(job.dir, '/mydir')
        self.assertEqual(job.meta.owner, 'bob')
        self.assertEqual(job.meta.description, 'My first description')
        self.assertEqual(job.meta.ms_filename, 'F123456.mzxml')
        self.assertEqual(job.meta.state, 'STOPPED')

        self.factory._makeJobDir.assert_called_with(job.id)
        self.factory._copyFile.assert_called_with(dbfile, job.id)
        o = mu.DBSession().query(mu.JobMeta.owner
                                 ).filter(mu.JobMeta.jobid == job.id).scalar()
        self.assertEqual(o, 'bob', 'job meta has been inserted')

    def test_fromid(self):
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        self.factory._makeJobSession = Mock(return_value=456)
        self.factory.id2jobdir = Mock(return_value=789)

        job = self.factory.fromId(jobid)

        self.assertEqual(job.owner, 'bob')
        self.assertEqual(job.db.session, 456)
        self.assertEqual(job.dir, 789)
        self.factory._makeJobSession.assert_called_with(jobid)
        self.factory.id2jobdir.assert_called_with(jobid)

    def test_fromid_notfoundindb(self):
        from magmaweb.job import JobNotFound
        jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        self.factory._makeJobSession = Mock()
        self.factory.id2jobdir = Mock()

        with self.assertRaises(JobNotFound) as exc:
            self.factory.fromId(jobid)
        self.assertEqual(exc.exception.jobid, jobid)
        self.assertEqual(exc.exception.message, "Job not found in database")

    def test_fromid_notfoundasdb(self):
        from magmaweb.job import JobNotFound
        jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        self.factory._getJobMeta = Mock(mu.JobMeta)

        with self.assertRaises(JobNotFound) as exc:
            self.factory.fromId(jobid)
        self.assertEqual(exc.exception.jobid, jobid)
        self.assertEqual(exc.exception.message, "Data of job not found")

    def test_submitQuery(self):
        self.factory.init_script = "# make magma available"
        job = self.factory.fromScratch('bob')

        self.factory.script_fn = 'script.sh'
        cmd = "magma add_structures -t smiles structures.dat results.db\n"
        jobquery = JobQuery(job.dir, cmd, ['structures.dat'])
        status_cb_url = 'http://example.com/status/{}.json'.format(job.id)
        jobquery.status_callback_url = status_cb_url

        self.factory.submitJob2Launcher = Mock()

        self.factory.submitQuery(jobquery, job)

        job_script_fn = os.path.join(jobquery.dir, self.factory.script_fn)
        job_script = open(job_script_fn).read()
        exp_script = "# make magma available\n"
        exp_script += "magma add_structures -t smiles "
        exp_script += "structures.dat results.db\n"
        self.assertMultiLineEqual(job_script, exp_script)
        jobmanager_query = {'jobdir': jobquery.dir + '/',
                            'executable': "/bin/sh",
                            'prestaged': [self.factory.script_fn,
                                          'results.db',
                                          'structures.dat'
                                          ],
                            "poststaged": ['results.db'],
                            "stderr": "stderr.txt",
                            "stdout": "stdout.txt",
                            "time_max": self.factory.time_max,
                            'arguments': [self.factory.script_fn],
                            'status_callback_url': status_cb_url
                            }
        self.factory.submitJob2Launcher.assert_called_with(jobmanager_query)
        self.assertEqual(job.state, 'INITIAL')

    def test_submitQuery_with_tarball(self):
        self.factory.tarball = 'Magma-1.1.tar.gz'
        self.factory.submitJob2Launcher = Mock()
        job = self.factory.fromScratch('bob')
        jobquery = JobQuery(job.dir, "", [])
        status_cb_url = 'http://example.com/status/{}.json'.format(job.id)
        jobquery.status_callback_url = status_cb_url

        self.factory.submitQuery(jobquery, job)

        jobmanager_query = {'jobdir': jobquery.dir + '/',
                            'executable': "/bin/sh",
                            'prestaged': [self.factory.script_fn,
                                          'results.db',
                                          # tarball is staged as well
                                          'Magma-1.1.tar.gz'
                                          ],
                            "poststaged": ['results.db'],
                            "stderr": "stderr.txt",
                            "stdout": "stdout.txt",
                            "time_max": self.factory.time_max,
                            'arguments': [self.factory.script_fn],
                            'status_callback_url': status_cb_url
                            }
        self.factory.submitJob2Launcher.assert_called_with(jobmanager_query)
        self.assertEqual(job.state, 'INITIAL')

    def test_submitQuery_no_jobmanager(self):
        from urllib2 import URLError
        from magmaweb.job import JobSubmissionError
        exc = URLError('[Errno 111] Connection refused')
        self.factory.submitJob2Launcher = Mock(side_effect=exc)
        job = self.factory.fromScratch('bob')
        jobquery = JobQuery(job.dir, "", [])
        status_cb_url = 'http://example.com/status/{}.json'.format(job.id)
        jobquery.status_callback_url = status_cb_url

        with self.assertRaises(JobSubmissionError):
            self.factory.submitQuery(jobquery, job)

        self.assertEqual(job.state, 'SUBMISSION_ERROR')

    @patch('urllib2.urlopen')
    def test_submitJob2Launcher(self, ua):
        import json
        body = {'foo': 'bar'}

        self.factory.submitJob2Launcher(body)

        # TODO replace with ua.assert_called_once_with(<hamcrest matcher>)
        req = ua.call_args[0][0]
        self.assertEquals(req.get_data(), json.dumps(body))
        self.assertEquals(req.get_full_url(), self.factory.submit_url)
        self.assertEquals(req.get_header('Content-type'), 'application/json')
        self.assertEquals(req.get_header('Accept'), 'application/json')

    def test_fromScratch(self):
        # mock/stub private methods which do external calls
        self.factory._makeJobDir = Mock(return_value='/mydir')
        self.factory._copyFile = Mock()
        db = initTestingDB(dataset=None)
        self.factory._makeJobSession = Mock(return_value=db)
        self.factory._addJobMeta = Mock()

        job = self.factory.fromScratch('bob')

        self.assertIsInstance(job.id, uuid.UUID)
        self.assertEqual(job.dir, '/mydir')
        self.assertEqual(job.meta.owner, 'bob')
        self.assertEqual(job.meta.state, 'STOPPED')
        self.assertEqual(job.meta.description, '')
        self.assertEqual(job.meta.ms_filename, '')
        self.assertEqual(job.db.maxMSLevel(), 0)
        self.factory._makeJobDir.assert_called_with(job.id)

    def test_cloneJob(self):
        oldjob = self.factory.fromScratch('ed')
        oldjob.description = 'My first description'
        oldjob.ms_filename = 'F123456.mzxml'

        job = self.factory.cloneJob(oldjob, 'bob')

        self.assertIsInstance(job.id, uuid.UUID)
        self.assertNotEqual(job.id, oldjob.id)
        self.assertEqual(job.meta.owner, 'bob')
        self.assertEqual(job.meta.description, 'My first description')
        self.assertEqual(job.meta.ms_filename, 'F123456.mzxml')
        self.assertEqual(job.meta.state, 'STOPPED')
        self.assertEqual(job.meta.parentjobid, oldjob.id)


class JobNotFound(unittest.TestCase):
    def test_it(self):
        from magmaweb.job import JobNotFound
        jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        e = JobNotFound('Job not found', jobid)
        self.assertEqual(e.jobid, jobid)
        self.assertEqual(e.message, 'Job not found')


class JobTestCase(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.parentjobid = uuid.UUID('22222222-2222-2222-2222-222222222222')
        self.jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        self.created_at = datetime.datetime(2012, 11, 14, 10, 48, 26, 504478)
        self.meta = mu.JobMeta(jobid=self.jobid,
                               description="My desc",
                               state='STOPPED',
                               parentjobid=self.parentjobid,
                               owner='bob',
                               created_at=self.created_at,
                               ms_filename='F1234.mzxml',
                               is_public=False,
                               )
        self.db = Mock(JobDb)
        self.jobdir = tempfile.mkdtemp()
        stderr = open(os.path.join(self.jobdir, 'stderr.txt'), 'w')
        stderr.write('Error log')
        stderr.close()
        self.job = Job(self.meta, self.jobdir, self.db)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.jobdir)

    def test_construct(self):
        self.assertEqual(self.job.db, self.db)
        self.assertEqual(self.job.dir, self.jobdir)

    def test_id(self):
        self.assertEqual(self.job.id, self.jobid)

    def test_name(self):
        self.assertEqual(self.job.__name__, str(self.jobid))

    def test_get_description(self):
        self.assertEqual(self.job.description, 'My desc')

    def test_set_description(self):
        self.job.description = 'My second description'

        self.assertEqual(self.job.description, 'My second description')
        self.assertEqual(self.job.meta.description, 'My second description')
        self.assertEqual(self.job.db.runInfo().description,
                         'My second description')

    def test_set_description_withoutruninfo(self):
        self.job.db.runInfo.return_value = None

        self.job.description = 'My second description'

        self.assertEqual(self.job.description, 'My second description')
        self.assertEqual(self.job.meta.description, 'My second description')

    def test_owner(self):
        self.assertEqual(self.job.owner, 'bob')

    def test_set_owner(self):
        self.job.owner = 'ed'
        self.assertEqual(self.job.owner, 'ed')

    def test_stderr(self):
        log = self.job.stderr()
        self.assertIsInstance(log, file)
        self.assertEqual(log.name, self.jobdir + '/stderr.txt')
        self.assertEqual(log.read(), 'Error log')

    def test_stderr_empty(self):
        open_stub = Mock(side_effect=IOError('File not found'))
        with patch('__builtin__.open', open_stub):
            log = self.job.stderr()
            self.assertEqual(log.read(), '')

    def test_jobquery(self):
        status_cb_url = 'http://example/status/{}.json'.format(self.job.id)
        jobquery = self.job.jobquery(status_cb_url, False)

        ejobquery = JobQuery(self.jobdir,
                             status_callback_url=status_cb_url,
                             restricted=False,
                             )

        self.assertIsInstance(jobquery, JobQuery)
        self.assertEqual(jobquery, ejobquery)

    def test_jobquery_restricted(self):
        status_cb_url = 'http://example/status/{}.json'.format(self.job.id)
        jobquery = self.job.jobquery(status_cb_url, True)

        ejobquery = JobQuery(self.jobdir,
                             status_callback_url=status_cb_url,
                             restricted=True,
                             )

        self.assertEqual(jobquery, ejobquery)

    def test_state(self):
        self.assertEquals(self.job.state, 'STOPPED')

    def test_set_state(self):
        self.job.state = 'RUNNING'

        self.assertEquals(self.job.state, 'RUNNING')

    def test_parent(self):
        self.assertEquals(self.job.parent, self.parentjobid)

    def test_set_parent(self):
        jid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        self.job.parent = jid

        self.assertEquals(self.job.parent, jid)

    def test_created_at(self):
        self.assertEqual(self.job.created_at, self.created_at)

    def test_ms_filename(self):
        self.assertEqual(self.job.ms_filename, 'F1234.mzxml')

    def test_set_ms_filename(self):
        self.job.ms_filename = 'F4567.mzxml'

        self.assertEqual(self.job.ms_filename, 'F4567.mzxml')
        self.assertEqual(self.job.meta.ms_filename, 'F4567.mzxml')
        self.assertEqual(self.job.db.runInfo().ms_filename, 'F4567.mzxml')

    def test_set_ms_filename_withoutruninfo(self):
        self.job.db.runInfo.return_value = None

        self.job.ms_filename = 'F4567.mzxml'

        self.assertEqual(self.job.ms_filename, 'F4567.mzxml')
        self.assertEqual(self.job.meta.ms_filename, 'F4567.mzxml')

    @patch('magmaweb.user.JobMeta')
    @patch('magmaweb.job.shutil')
    def test_delete(self, shutl, jm):
        self.job.db.session = Mock()

        self.job.delete()

        self.job.db.session.remove.assert_called_with()
        shutl.rmtree.assert_called_with(self.jobdir)
        jm.delete.assert_called_with(self.meta)

    def test_is_public(self):
        self.assertEquals(self.job.is_public, False)

    def test_set_is_public(self):
        self.job.is_public = True

        self.assertEquals(self.job.is_public, True)


class JobDbTestCaseAbstract(unittest.TestCase):
    def setUp(self):
        # mock job session
        self.session = initTestingDB()
        self.job = JobDb(self.session)


class JobDbTestCase(JobDbTestCaseAbstract):
    def test_construct(self):
        self.assertEqual(self.job.session, self.session)

    def test_runInfo(self):
        runInfo = self.job.runInfo()
        self.assertEqual(runInfo.n_reaction_steps, 2)
        self.assertEqual(runInfo.metabolism_types, "phase1,phase2")

        self.assertEqual(runInfo.ms_filename, 'F123456.mzxml')
        self.assertEqual(runInfo.abs_peak_cutoff, 1000)
        self.assertEqual(runInfo.max_ms_level, 3)
        self.assertEqual(runInfo.precursor_mz_precision, 10)

        self.assertEqual(runInfo.ionisation_mode, -1)
        self.assertEqual(runInfo.max_broken_bonds, 4)
        self.assertEqual(runInfo.ms_intensity_cutoff, 200000.0)
        self.assertEqual(runInfo.msms_intensity_cutoff, 50)
        self.assertEqual(runInfo.mz_precision, 10)
        self.assertEqual(runInfo.mz_precision_abs, 0.002)
        self.assertEqual(runInfo.description, 'My first description')
        self.assertEqual(runInfo.max_water_losses, 1)

    def test_runInfo_maxrunid(self):
        self.session.add(Run(
            n_reaction_steps=2, metabolism_types='phase1,phase2',
            ionisation_mode=-1, skip_fragmentation=True,
            ms_intensity_cutoff=200000.0, msms_intensity_cutoff=50,
            mz_precision=10, mz_precision_abs=0.002, use_all_peaks=True,
            ms_filename='F123456.mzxml', abs_peak_cutoff=1000,
            max_ms_level=3, precursor_mz_precision=10,
            max_broken_bonds=4, description='My second description',
            max_water_losses=1,
        ))

        runInfo = self.job.runInfo()

        # run with highest id is returned
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
        expected_chromatogram = {'scans': [{'id': 641, 'rt': 933.317,
                                            'intensity': 807577.0, 'ap': 0},
                                           {'id': 870, 'rt': 1254.15,
                                            'intensity': 1972180.0, 'ap': 0},
                                           ],
                                 'cutoff': 200000.0,
                                 }
        self.assertEqual(self.job.chromatogram(), expected_chromatogram)

    def test_chromatogram_with_assigned_peaks(self):
        metid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz, metid)

        expected_chromatogram = {'scans': [{'id': 641, 'rt': 933.317,
                                            'intensity': 807577.0, 'ap': 1},
                                           {'id': 870, 'rt': 1254.15,
                                            'intensity': 1972180.0, 'ap': 0},
                                           ],
                                 'cutoff': 200000.0,
                                 }
        self.assertEqual(self.job.chromatogram(), expected_chromatogram)

    def test_metabolitesTotalCount(self):
        self.assertEqual(self.job.metabolitesTotalCount(), 2)

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


class JobDbEmptyDatasetTestCase(unittest.TestCase):
    def setUp(self):
        # mock job session
        self.session = initTestingDB(dataset=None)
        self.job = JobDb(self.session)

    def test_runInfo(self):
        runInfo = self.job.runInfo()
        self.assertIsNone(runInfo)

    def test_maxMSLevel(self):
        maxmslevel = self.job.maxMSLevel()
        self.assertEqual(maxmslevel, 0)

    def test_chromatogram(self):
        expected_chromatogram = {'scans': [], 'cutoff': None}
        self.assertEqual(self.job.chromatogram(), expected_chromatogram)

    def test_metabolitesTotalCount(self):
        self.assertEqual(self.job.metabolitesTotalCount(), 0)


class JobDbMetabolitesTestCase(JobDbTestCaseAbstract):
    def test_default(self):
        response = self.job.metabolites()
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += '?cid=152432">CID: 152432</a>'
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
                    'reactionsequence': [u'PARENT'],
                    'smiles': u'Oc1ccccc1O',
                    'mim': 110.03677, 'logp':1.231,
                    'assigned': False,
                    'reference': url1
                }, {
                    'isquery': True, 'level': 0, 'metid': 352,
                    'mol': u"Molfile of dihydroxyphenyl-valerolactone",
                    'molformula': u"C11H12O4",
                    'nhits': None,
                    'nhits': 1,
                    'origin': u"dihydroxyphenyl-valerolactone",
                    'probability': 1,
                    'reactionsequence': [u"PARENT", u"CHILD"],
                    'smiles': u"O=C1OC(Cc2ccc(O)c(O)c2)CC1",
                    'mim': 208.07355, 'logp':2.763,
                    'assigned': False,
                    'reference': url2
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


class JobDbMetabolites2csvTestCase(JobDbTestCaseAbstract):
    def test_it(self):
        csvfile = self.job.metabolites2csv(self.job.metabolites()['rows'])
        import csv
        import StringIO
        expected_csvfile = StringIO.StringIO()
        cols = ['origin', 'smiles', 'probability', 'reactionsequence',
                'nhits', 'molformula', 'mim', 'isquery', 'logp', 'reference']
        csvwriter = csv.DictWriter(expected_csvfile, cols)
        csvwriter.writeheader()
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += '?cid=152432">CID: 152432</a>'
        csvwriter.writerow({'origin': 'pyrocatechol',
                            'smiles': 'Oc1ccccc1O',
                            'probability': 1.0,
                            'reactionsequence': 'PARENT',
                            'nhits': 1,
                            'molformula': 'C6H6O2',
                            'isquery': True,
                            'mim': 110.03677,
                            'logp': 1.231,
                            'reference': url1,
                            })
        csvwriter.writerow({'origin': 'dihydroxyphenyl-valerolactone',
                            'smiles': 'O=C1OC(Cc2ccc(O)c(O)c2)CC1',
                            'probability': 1.0,
                            'reactionsequence': 'PARENT|CHILD',
                            'nhits': 1,
                            'molformula': 'C11H12O4',
                            'isquery': True,
                            'mim': 208.07355,
                            'logp': 2.763,
                            'reference': url2,
                            })
        self.assertMultiLineEqual(csvfile.getvalue(),
                                  expected_csvfile.getvalue())

    def test_with_score(self):
        mets = self.job.metabolites(scanid=641)['rows']
        csvfile = self.job.metabolites2csv(mets)
        import csv
        import StringIO
        expected_csvfile = StringIO.StringIO()
        cols = ['origin', 'smiles', 'probability', 'reactionsequence',
                'nhits', 'molformula', 'mim', 'isquery', 'logp', 'reference',
                'score']
        csvwriter = csv.DictWriter(expected_csvfile, cols)
        csvwriter.writeheader()
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        csvwriter.writerow({'origin': 'pyrocatechol', 'smiles': 'Oc1ccccc1O',
                            'probability': 1.0, 'reactionsequence': 'PARENT',
                            'nhits': 1, 'molformula': 'C6H6O2',
                            'isquery': True, 'score': 200.0, 'mim': 110.03677,
                            'logp': 1.231,
                            'reference': url1
                            })
        self.assertMultiLineEqual(csvfile.getvalue(),
                                  expected_csvfile.getvalue())

    def test_some_columns(self):
        cols = ['origin', 'mim']
        mets = self.job.metabolites(scanid=641)['rows']

        response = self.job.metabolites2csv(mets, cols=cols)

        self.assertEquals(
            response.getvalue(),
            'origin,mim\r\n' +
            'pyrocatechol,110.03677\r\n'
        )


class JobMetabolites2sdfTestCase(JobDbTestCaseAbstract):
    def test_it(self):
        sdffile = self.job.metabolites2sdf(self.job.metabolites()['rows'])
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += '?cid=152432">CID: 152432</a>'

        expected_sdf = """Molfile> <origin>
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
{url1}

$$$$
Molfile of dihydroxyphenyl-valerolactone> <origin>
dihydroxyphenyl-valerolactone

> <smiles>
O=C1OC(Cc2ccc(O)c(O)c2)CC1

> <probability>
1.0

> <reactionsequence>
PARENT
CHILD

> <nhits>
1

> <molformula>
C11H12O4

> <mim>
208.07355

> <logp>
2.763

> <reference>
{url2}

$$$$
""".format(url1=url1, url2=url2)

        self.assertMultiLineEqual(sdffile, expected_sdf)

    def test_with_scan(self):
        """Include score prop"""
        mets = self.job.metabolites(scanid=641)['rows']
        sdffile = self.job.metabolites2sdf(mets)
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        expected_sdf = """Molfile> <origin>
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
{url1}

> <score>
200.0

$$$$
""".format(url1=url1)

        self.assertMultiLineEqual(sdffile, expected_sdf)

    def test_some_columns(self):
        cols = ['origin', 'mim']

        mets = self.job.metabolites(scanid=641)['rows']
        sdffile = self.job.metabolites2sdf(mets, cols=cols)

        expected_sdf = """Molfile> <origin>
pyrocatechol

> <mim>
110.03677

$$$$
"""

        self.assertMultiLineEqual(sdffile, expected_sdf)


class JobScansWithMetabolitesTestCase(JobDbTestCaseAbstract):
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

    def test_filteron_deltappm(self):
        filters = [{"type": "numeric", "comparison": "gt",
                    "value": 0, "field": "nhits"},
                   {"type": "numeric", "value": "200",
                    "comparison": "lt", "field": "deltappm"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])


class JobMSpectraTestCase(JobDbTestCaseAbstract):
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
            'cutoff': 13950500.0,
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


class JobFragmentsTestCase(JobDbTestCaseAbstract):
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
                'formula': u'C5H4',
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
                    'formula': 'C3H5O3',
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
                    'score': 65,
                    'formula': 'C3H5O3',
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
                'isAssigned': False,
                'formula': 'C3H5O3',
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
            'deltappm': 5.0781684060061766e-08,
            'formula': 'C3H5O3',
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
                'formula': "C5H4",
            }], 'expanded': True
        })

    def test_badfragment(self):
        from magmaweb.job import FragmentNotFound
        with self.assertRaises(FragmentNotFound):
            self.job.fragments(metid=70002, scanid=641, node='root')


class JobWithAllPeaksTestCase(unittest.TestCase):
    def setUp(self):
        self.job = JobDb(initTestingDB(dataset='useallpeaks'))

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
                    'reactionsequence': [u'sulfation_(aromatic_hydroxyl)'],
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
                'isAssigned': False,
                'formula': "C4H6O2",
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
                'formula': "C4H6O2",
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
                    'deltappm': 2.3630267823129625e-12,
                    'formula': "C4H6O2",
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
                    'deltappm': -7.021641217476223e-09,
                    'formula': "C4H6O2",
                }]
            }], 'expanded': True
        })
