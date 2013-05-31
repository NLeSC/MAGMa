"""Tests for magmaweb.job except JobQuery and JobDb"""
import datetime
import os
import uuid
import unittest
from mock import Mock, patch
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import transaction
from magmaweb.job import JobFactory, Job, JobDb, make_job_factory, JobQuery
from magmaweb.job import JobError, JobIncomplete, MissingDataError
from magmaweb.job import JobNotFound
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
    url1 = u'<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
    url1 += u'?cid=289">CID: 289</a>'
    url2 = u'<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
    url2 += u'?cid=152432">CID: 152432</a>'
    session.add(Run(
        n_reaction_steps=2, metabolism_types=u'phase1,phase2',
        ionisation_mode=-1, skip_fragmentation=True,
        ms_intensity_cutoff=200000.0, msms_intensity_cutoff=50,
        mz_precision=10, mz_precision_abs=0.002, use_all_peaks=True,
        ms_filename=u'F123456.mzxml', abs_peak_cutoff=1000,
        max_ms_level=3, precursor_mz_precision=10,
        max_broken_bonds=4, description=u'My first description',
        max_water_losses=1,
    ))
    session.add(Metabolite(
        metid=72, mol=u'Molfile', level=0, probability=1.0,
        reactionsequence=['PARENT'], smiles=u'Oc1ccccc1O',
        molformula=u'C6H6O2', isquery=True, nhits=1,
        origin=u'pyrocatechol', mim=110.03677, logp=1.231,
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
        atoms=u"0,1,2,3,4,5,6,7",
        deltah=-1.0,
        deltappm=-1.84815979523607e-08,
        formula=u"C5H4",
    ))
    # fragments of metid=352 + scanid=870
    session.add(Metabolite(
        isquery=True, level=0, metid=352,
        mol=u"Molfile of dihydroxyphenyl-valerolactone",
        molformula=u"C11H12O4",
        origin=u"dihydroxyphenyl-valerolactone",
        probability=1.0, reactionsequence=u"PARENT\nCHILD\n",
        smiles=u"O=C1OC(Cc2ccc(O)c(O)c2)CC1",
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
        atoms=u"0,1,2,3,4,5,6,7,8,9,10,11,12,13,14", deltah=-1,
        deltappm=-8.18675317722029e-09,
        formula=u"C3H5O3",
    ), Fragment(
        fragid=1708, metid=352, scanid=871, mass=123.0446044689,
        mz=123.04508972167969, score=201, parentfragid=1707,
        atoms=u"6,7,8,9,10,11,12,13,14", deltah=0,
        deltappm=3.943698856902144e-12,
        formula=u"C3H5O3",
    ), Fragment(
        fragid=1709, metid=352, scanid=871, mass=164.08372962939995,
        mz=163.07623291015625, score=65, parentfragid=1707,
        atoms=u"3,4,5,6,7,8,9,10,11,12,13,14", deltah=-1,
        deltappm=-1.235815738001507e-08,
        formula=u"C3H5O3",
    ), Fragment(
        fragid=1710, scanid=872, metid=352, mass=116.0626002568,
        mz=119.08654022216797, score=4, parentfragid=1709,
        atoms=u"4,5,6,7,8,9,11,13,14", deltah=3,
        deltappm=5.0781684060061766e-08,
        formula=u"C3H5O3",
    )])

    session.flush()


def populateWithUseAllPeaks(session):
    """ Dataset with multiple fragments of same metabolite on lvl1 scan """
    session.add(Run(
        n_reaction_steps=2, metabolism_types=u'phase1,phase2',
        ionisation_mode=-1, skip_fragmentation=True,
        ms_intensity_cutoff=200000.0, msms_intensity_cutoff=50,
        mz_precision=10, mz_precision_abs=0.002, use_all_peaks=False,
        ms_filename=u'F123456.mzxml', abs_peak_cutoff=1000,
        max_ms_level=3, precursor_mz_precision=10,
        max_broken_bonds=4, description=u'My second description',
        max_water_losses=1,
    ))
    session.add(Metabolite(
        metid=12,
        level=1,
        probability=0.119004,
        reactionsequence=['sulfation_(aromatic_hydroxyl)'],
        smiles=u'Oc1ccc(CC2OC(=O)CC2)cc1OS(O)(=O)=O',
        molformula=u'C11H12O7S',
        isquery=False,
        origin=u'5-(3,4)-dihydroxyphenyl-g-valerolactone (F)',
        mol=u'Molfile',
        reference=u'',
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
        atoms=u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
        deltah=-1.0,
        deltappm=-7.046696487857745e-09,
        formula=u"C4H6O2",
    ), Fragment(
        fragid=18,
        metid=12,
        scanid=1,
        mz=287.023132324219,
        mass=288.0303734299,
        score=0.5,
        parentfragid=0,
        atoms=u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
        deltah=-1.0,
        deltappm=-7.020570507205176e-09,
        formula=u"C4H6O2",
    ), Fragment(
        fragid=19,
        metid=12,
        scanid=2,
        mz=207.066223144531,
        mass=207.0657338415,
        score=0.5,
        parentfragid=18,
        atoms=u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14',
        deltah=0.0,
        deltappm=2.3630267823129625e-12,
        formula=u"C4H6O2",
    ), Fragment(
        fragid=20,
        metid=12,
        scanid=2,
        mz=287.022827148438,
        mass=288.0303734299,
        score=0.0,
        parentfragid=18,
        atoms=u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
        deltah=-1.0,
        deltappm=-7.021641217476223e-09,
        formula=u"C4H6O2",
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
        self.factory = JobFactory(root_dir=self.root_dir,
                                  submit_url='http://localhost:9998/job',
                                  )

        # fill user db
        transaction.begin()
        engine = create_engine('sqlite:///:memory:')
        mu.DBSession.configure(bind=engine)
        mu.Base.metadata.create_all(engine)  # @UndefinedVariable
        mu.DBSession().add(mu.User(u'bob', u'Bob', u'bob@example.com'))
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        mu.DBSession().add(mu.JobMeta(jobid, owner=u'bob'))

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

        job = self.factory.fromDb(dbfile, u'bob')

        self.assertIsInstance(job.id, uuid.UUID)
        self.assertEqual(job.dir, u'/mydir')
        self.assertEqual(job.meta.owner, u'bob')
        self.assertEqual(job.meta.description, u'My first description')
        self.assertEqual(job.meta.ms_filename, u'F123456.mzxml')
        self.assertEqual(job.meta.state, u'STOPPED')

        self.factory._makeJobDir.assert_called_with(job.id)
        self.factory._copyFile.assert_called_with(dbfile, job.id)
        o = mu.DBSession().query(mu.JobMeta.owner
                                 ).filter(mu.JobMeta.jobid == job.id).scalar()
        self.assertEqual(o, u'bob', 'job meta has been inserted')

    def test_fromid(self):
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        self.factory._makeJobSession = Mock(return_value=456)
        self.factory.id2jobdir = Mock(return_value=789)

        job = self.factory.fromId(jobid)

        self.assertEqual(job.owner, u'bob')
        self.assertEqual(job.db.session, 456)
        self.assertEqual(job.dir, 789)
        self.factory._makeJobSession.assert_called_with(jobid)
        self.factory.id2jobdir.assert_called_with(jobid)

    def test_fromid_notfoundindb(self):
        jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        self.factory._makeJobSession = Mock()
        self.factory.id2jobdir = Mock()

        with self.assertRaises(JobNotFound) as exc:
            self.factory.fromId(jobid)
        self.assertEqual(exc.exception.jobid, jobid)
        self.assertEqual(exc.exception.message, "Job not found in database")

    def test_fromid_notfoundasdb(self):
        jobid = uuid.UUID('11111111-1111-1111-1111-111111111111')
        self.factory._getJobMeta = Mock(mu.JobMeta)

        with self.assertRaises(JobNotFound) as exc:
            self.factory.fromId(jobid)
        self.assertEqual(exc.exception.jobid, jobid)
        self.assertEqual(exc.exception.message, "Data of job not found")

    def test_submitQuery(self):
        self.factory.init_script = "# make magma available"
        job = self.factory.fromScratch(u'bob')

        self.factory.script_fn = 'script.sh'
        cmd = "magma add_structures -t smiles structures.dat results.db\n"
        jobquery = JobQuery(job.dir, cmd, ['structures.dat'])
        status_cb_url = 'http://example.com/status/{}.json'.format(job.id)
        jobquery.status_callback_url = status_cb_url

        launcher_url = 'http://localhost:9998/job/'
        launcher_url += '70a00fe2-f698-41ed-b28c-b37c22f10440'
        self.factory.submitJob2Launcher = Mock(return_value=launcher_url)

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
                            'arguments': [self.factory.script_fn],
                            'status_callback_url': status_cb_url
                            }
        self.factory.submitJob2Launcher.assert_called_with(jobmanager_query)
        self.assertEqual(job.state, u'INITIAL')
        self.assertEqual(job.launcher_url, launcher_url)

    def test_submitQuery_with_tarball(self):
        self.factory.tarball = 'Magma-1.1.tar.gz'
        self.factory.submitJob2Launcher = Mock()
        job = self.factory.fromScratch(u'bob')
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
                            'arguments': [self.factory.script_fn],
                            'status_callback_url': status_cb_url
                            }
        self.factory.submitJob2Launcher.assert_called_with(jobmanager_query)
        self.assertEqual(job.state, u'INITIAL')

    def test_submitQuery_no_joblauncher(self):
        from requests.exceptions import ConnectionError
        from magmaweb.job import JobSubmissionError
        exc = ConnectionError('[Errno 111] Connection refused')
        self.factory.submitJob2Launcher = Mock(side_effect=exc)
        job = self.factory.fromScratch(u'bob')
        jobquery = JobQuery(job.dir, "", [])
        status_cb_url = 'http://example.com/status/{}.json'.format(job.id)
        jobquery.status_callback_url = status_cb_url

        with self.assertRaises(JobSubmissionError):
            self.factory.submitQuery(jobquery, job)

        self.assertEqual(job.state, 'SUBMISSION_ERROR')

    def test_submitQuery_invalid_joblauncher(self):
        from requests.exceptions import HTTPError
        from magmaweb.job import JobSubmissionError
        exc = HTTPError('422 Client Error: Unprocessable Entity')
        self.factory.submitJob2Launcher = Mock(side_effect=exc)
        job = self.factory.fromScratch(u'bob')
        jobquery = JobQuery(job.dir, "", [])
        status_cb_url = 'http://example.com/status/{}.json'.format(job.id)
        jobquery.status_callback_url = status_cb_url

        with self.assertRaises(JobSubmissionError):
            self.factory.submitQuery(jobquery, job)

        self.assertEqual(job.state, 'SUBMISSION_ERROR')

    @patch('requests.post')
    def test_submitJob2Launcher(self, ua):
        from requests import Response
        create_url = 'http://localhost:9998/job/'
        create_url += '11111111-1111-1111-1111-111111111111'
        uresp = Response()
        uresp.status_code = 204
        uresp.headers = {'Location': create_url}
        ua.return_value = uresp
        body = {'foo': 'bar'}

        resp = self.factory.submitJob2Launcher(body)

        headers = {'Content-Type': 'application/json',
                   'Accept': 'application/json',
                   }
        ua.assert_called_with('http://localhost:9998/job',
                              data='{"foo": "bar"}',
                              headers=headers,
                              )

        launcher_url = 'http://localhost:9998/job/'
        launcher_url += '11111111-1111-1111-1111-111111111111'
        self.assertEquals(resp, launcher_url)

    def test_fromScratch(self):
        # mock/stub private methods which do external calls
        self.factory._makeJobDir = Mock(return_value='/mydir')
        self.factory._copyFile = Mock()
        db = initTestingDB(dataset=None)
        self.factory._makeJobSession = Mock(return_value=db)
        self.factory._addJobMeta = Mock()

        job = self.factory.fromScratch(u'bob')

        self.assertIsInstance(job.id, uuid.UUID)
        self.assertEqual(job.dir, u'/mydir')
        self.assertEqual(job.meta.owner, u'bob')
        self.assertEqual(job.meta.state, u'STOPPED')
        self.assertEqual(job.meta.description, u'')
        self.assertEqual(job.meta.ms_filename, u'')
        self.assertEqual(job.db.maxMSLevel(), 0)
        self.factory._makeJobDir.assert_called_with(job.id)

    def test_cloneJob(self):
        oldjob = self.factory.fromScratch(u'ed')
        oldjob.description = u'My first description'
        oldjob.ms_filename = u'F123456.mzxml'

        job = self.factory.cloneJob(oldjob, u'bob')

        self.assertIsInstance(job.id, uuid.UUID)
        self.assertNotEqual(job.id, oldjob.id)
        self.assertEqual(job.meta.owner, u'bob')
        self.assertEqual(job.meta.description, u'My first description')
        self.assertEqual(job.meta.ms_filename, u'F123456.mzxml')
        self.assertEqual(job.meta.state, u'STOPPED')
        self.assertEqual(job.meta.parentjobid, oldjob.id)

    @patch('requests.delete')
    def test_cancel(self, ua):
        job = self.factory.fromScratch(u'ed')
        url = 'http://localhost:9998/job/70a00fe2-f698-41ed-b28c-b37c22f10440'
        job.launcher_url = url

        self.factory.cancel(job)

        ua.assert_called_with(url)


class JobNotFoundTestCase(unittest.TestCase):
    def test_it(self):
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
                               description=u"My desc",
                               state=u'STOPPED',
                               parentjobid=self.parentjobid,
                               owner=u'bob',
                               created_at=self.created_at,
                               ms_filename=u'F1234.mzxml',
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
        self.assertEqual(self.job.description, u'My desc')

    def test_set_description(self):
        self.job.description = u'My second description'

        self.assertEqual(self.job.description, u'My second description')
        self.assertEqual(self.job.meta.description, u'My second description')
        self.assertEqual(self.job.db.runInfo().description,
                         u'My second description')

    def test_set_description_withoutruninfo(self):
        self.job.db.runInfo.return_value = None

        self.job.description = u'My second description'

        self.assertEqual(self.job.description, u'My second description')
        self.assertEqual(self.job.meta.description, u'My second description')

    def test_owner(self):
        self.assertEqual(self.job.owner, u'bob')

    def test_set_owner(self):
        self.job.owner = u'ed'
        self.assertEqual(self.job.owner, u'ed')

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
        self.assertEquals(self.job.state, u'STOPPED')

    def test_set_state(self):
        self.job.state = u'RUNNING'

        self.assertEquals(self.job.state, u'RUNNING')

    def test_parent(self):
        self.assertEquals(self.job.parent, self.parentjobid)

    def test_set_parent(self):
        jid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        self.job.parent = jid

        self.assertEquals(self.job.parent, jid)

    def test_created_at(self):
        self.assertEqual(self.job.created_at, self.created_at)

    def test_ms_filename(self):
        self.assertEqual(self.job.ms_filename, u'F1234.mzxml')

    def test_set_ms_filename(self):
        self.job.ms_filename = u'F4567.mzxml'

        self.assertEqual(self.job.ms_filename, u'F4567.mzxml')
        self.assertEqual(self.job.meta.ms_filename, u'F4567.mzxml')
        self.assertEqual(self.job.db.runInfo().ms_filename, u'F4567.mzxml')

    def test_set_ms_filename_withoutruninfo(self):
        self.job.db.runInfo.return_value = None

        self.job.ms_filename = u'F4567.mzxml'

        self.assertEqual(self.job.ms_filename, u'F4567.mzxml')
        self.assertEqual(self.job.meta.ms_filename, u'F4567.mzxml')

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

    def test_is_complete(self):
        self.job.state = u'STOPPED'
        self.assertTrue(self.job.is_complete())

    def test_is_complete_running(self):
        running_states = (u'INITIAL', u'PRE_STAGING', u'RUNNING',
                          u'POST_STAGING', u'Progress: 50%')
        for running_state in running_states:
            self.job.state = running_state
            with self.assertRaises(JobIncomplete) as e:
                self.job.is_complete()

            self.assertEqual(e.exception.job, self.job)

    def test_is_complete_error(self):
        self.job.state = u'ERROR'
        with self.assertRaises(JobError) as e:
            self.job.is_complete()

        self.assertEqual(e.exception.job, self.job)

    def test_is_complete_filled(self):
        self.job.db.hasMolecules.return_value = True
        self.job.db.hasMspectras.return_value = True
        self.job.db.hasFragments.return_value = True

        complete = self.job.is_complete(True)

        self.assertTrue(complete)

    def test_is_complete_unfilled_molecules(self):
        self.job.db.hasMolecules = Mock(return_value=False)
        self.job.db.hasMspectras = Mock(return_value=True)
        self.job.db.hasFragments = Mock(return_value=True)

        with self.assertRaises(MissingDataError) as e:
            self.job.is_complete(True)

        self.assertEqual(e.exception.job, self.job)
        self.assertEqual(e.exception.message, 'No molecules found')

    def test_is_complete_unfilled_mspectra(self):
        self.job.db.hasMolecules = Mock(return_value=True)
        self.job.db.hasMspectras = Mock(return_value=False)
        self.job.db.hasFragments = Mock(return_value=True)

        with self.assertRaises(MissingDataError) as e:
            self.job.is_complete(True)

        self.assertEqual(e.exception.job, self.job)
        self.assertEqual(e.exception.message, 'No mass spectras found')

    def test_is_complete_unfilled_fragments(self):
        self.job.db.hasMolecules = Mock(return_value=True)
        self.job.db.hasMspectras = Mock(return_value=True)
        self.job.db.hasFragments = Mock(return_value=False)

        with self.assertRaises(MissingDataError) as e:
            self.job.is_complete(True)

        self.assertEqual(e.exception.job, self.job)
        self.assertEqual(e.exception.message, 'No fragments found')

    def test_launcher_url(self):
        url = 'http://localhost:9998/job/70a00fe2-f698-41ed-b28c-b37c22f10440'
        self.job.launcher_url = url

        self.assertEqual(self.job.launcher_url, url)
        self.assertEqual(self.job.meta.launcher_url, url)
