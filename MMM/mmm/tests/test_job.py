import unittest
from mmm.job import JobFactory, Job
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

def initTestingDB(url = 'sqlite://'):
    """Creates testing db and populates with test data"""
    engine = create_engine(url) # default in memory db
    session = sessionmaker(bind=engine)
    dbh = session()
    from mmm.models import Base
    Base.metadata.create_all(engine)
    populateTestingDB(dbh)
    return dbh

def populateTestingDB(session):
    """Polulates test db with data

    Adds 1 metabolite with one fragment.
    Adds 1 metabolite with one fragment which has 2 child fragments of which one has another child fragment.
    Run, Scan and Peak are filled to create a working run.

    session
        session connection to db
    """
    from mmm.models import Metabolite, Scan, Peak, Fragment, Run
    session.add(Run(
        n_reaction_steps=2, use_phase1=True, use_phase2=True,
        ionisation=-1, use_fragmentation=True,
        ms_intensity_cutoff=200000.0, msms_intensity_cutoff=0.5,
        mz_precision=0.01, use_msms_only=True
    ))
    session.add(Metabolite(
        metid=72, mol='Molfile', level=0, probability=1.0,
        reactionsequence='PARENT', smiles='Oc1ccccc1O',
        molformula='C6H6O2', isquery=True,
        origin='pyrocatechol'
    ))
    session.add(Scan(
        scanid=641, mslevel=1, rt=933.317, lowmz=90.3916, highmz=1197.78,
        basepeakmz=305.034, basepeakintensity=807577.0, totioncurrent=5957620
    ))
    session.add_all([
        Peak(scanid=641, mz=305.033508300781, intensity=807576.625), # basepeak
        Peak(scanid=641, mz=109.0295639038086, intensity=345608.65625) # peak of metabolite
    ])
    session.add(Fragment(
        fragid=948,
        metid=72,
        scanid=641,
        mz=109.0296783,
        mass=110.0367794368,
        score=200,
        parentfragid=0,
        atoms="0,1,2,3,4,5,6,7",
        deltah=-1.0
    ))
    # fragments of metid=352 + scanid=870
    session.add(Metabolite(
        isquery=True, level=0, metid=352, mol="Molfile of dihydroxyphenyl-valerolactone",
        molformula="C11H12O4",
        origin="dihydroxyphenyl-valerolactone",
        probability=1.0, reactionsequence="PARENT",
        smiles="O=C1OC(Cc2ccc(O)c(O)c2)CC1"
    ))
    session.add_all([Scan(
        scanid=870, mslevel=1, rt=1254.15, lowmz=91.0302, highmz=1171.51,
        basepeakmz=287.023, basepeakintensity=1972180.0, totioncurrent=9265290
    ), Scan(
        scanid=871, mslevel=2, rt=1254.93, lowmz=51.5211, highmz=216.864,
        basepeakmz=163.076, basepeakintensity=279010.0, totioncurrent=809307,
        precursormz=207.0663147, precursorintensity=293096.0, precursorscanid=870
    ), Scan(
        scanid=872, mslevel=3, rt=1256.77, lowmz=50.3338, highmz=172.155,
        basepeakmz=119.087, basepeakintensity=17387.0, totioncurrent=236842,
        precursormz=163.0762329, precursorintensity=6163.73, precursorscanid=871
    )])
    session.add_all([
        Peak(scanid=870, mz=287.022979736328, intensity=1972180.625), # basepeak
        Peak(scanid=870, mz=207.066284179688, intensity=293095.84375), # peak of metabolite
        Peak(scanid=871, mz=163.076232910156, intensity=279010.28125), # basepeak and peak of frag 1709
        Peak(scanid=871, mz=123.04508972168, intensity=211603.046875), # peak of frag 1708
        Peak(scanid=872, mz=119.086540222168, intensity=17386.958984375), # basepeak and peak of frag 1710
    ])
    session.add_all([Fragment(
        fragid=1707, metid=352, scanid=870, mass=208.0735588736,
        mz=207.0663147, score=100, parentfragid=0,
        atoms="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14", deltah=-1
    ),Fragment(
        fragid=1708, metid=352, scanid=871, mass=123.0446044689,
        mz=123.04508972167969, score=201, parentfragid=1707,
        atoms="6,7,8,9,10,11,12,13,14", deltah=0
    ),Fragment(
        fragid=1709, metid=352, scanid=871, mass=164.08372962939995,
        mz=163.07623291015625, score=65, parentfragid=1707,
        atoms="3,4,5,6,7,8,9,10,11,12,13,14", deltah=-1
    ),Fragment(
        fragid=1710, scanid=872, metid=352, mass=116.0626002568,
        mz=119.08654022216797, score=4, parentfragid=1709,
        atoms="4,5,6,7,8,9,11,13,14", deltah=3
    )])
    session.flush()

class JobFactoryTestCase(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.jobrootdir = tempfile.mkdtemp()
        self.dbname = 'results.db'
        self.factory = JobFactory(self.jobrootdir, self.dbname)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.jobrootdir)

    def test_hasrootdir(self):
        self.assertEqual(self.factory.jobrootdir, self.jobrootdir)

    def test_hasdbname(self):
        self.assertEqual(self.factory.dbname, self.dbname)

    def test_id2url(self):
        import uuid, os
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jobdbfn = 'sqlite:///'+os.path.join(self.jobrootdir, str(jobid), self.dbname)
        self.assertEqual(self.factory.id2url(jobid), jobdbfn)

    def test_id2jobdir(self):
        import uuid, os
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jobdir = os.path.join(self.jobrootdir, str(jobid))
        self.assertEqual(self.factory.id2jobdir(jobid), jobdir)

    def test_id2db(self):
        import uuid, os
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jobdbfn = os.path.join(self.jobrootdir, str(jobid), self.dbname)
        self.assertEqual(self.factory.id2db(jobid), jobdbfn)

    def test_fromquery(self):
        import os, uuid
        dbfile = os.tmpfile()
        dbfile.write('bla')

        job = self.factory.fromQuery(dbfile)

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
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        job = self.factory.fromId(jobid)
        self.assertIsInstance(job, Job)
        self.assertEqual(job.id, jobid)
        self.assertEqual(
            str(job.session.get_bind().url),
            self.factory.id2url(jobid),
            'job has dbsession set to sqlite file in job dir'
        )

class JobTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.jobid = uuid.uuid1()
        # mock job session
        self.session = initTestingDB()
        self.job = Job(self.jobid, self.session)

    def test_construct(self):
        self.assertEqual(self.job.id, self.jobid)
        self.assertEqual(self.job.session, self.session)

    def test_runInfo(self):
        runInfo = self.job.runInfo()
        self.assertEqual(runInfo.n_reaction_steps, 2)
        self.assertEqual(runInfo.use_phase1, True)
        self.assertEqual(runInfo.use_phase2, True)
        self.assertEqual(runInfo.ionisation, u'-1')
        self.assertEqual(runInfo.use_fragmentation, True)
        self.assertEqual(runInfo.ms_intensity_cutoff, 200000.0)
        self.assertEqual(runInfo.msms_intensity_cutoff, 0.5)
        self.assertEqual(runInfo.mz_precision, 0.01)
        self.assertEqual(runInfo.use_msms_only, True)

    def test_maxMSLevel(self):
        maxmslevel = self.job.maxMSLevel()
        self.assertEqual(maxmslevel, 3)

    def test_extractedionchromatogram(self):
        metid = 72
        eic = self.job.extractedIonChromatogram(metid)
        self.assertEqual(eic,[{
            'rt': 933.317,
            'intensity': 345608.65625
        },{
            'rt': 1254.15,
            'intensity': 0
        }])

    def test_chromatogram(self):
        self.assertEqual(
            self.job.chromatogram(),[
            { 'id': 641, 'rt': 933.317, 'intensity': 807577.0 },
            { 'id': 870, 'rt': 1254.15, 'intensity': 1972180.0 }
        ])

class JobMetabolitesTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB())

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
                    'nr_scans': 1,
                    'origin': u'pyrocatechol',
                    'probability': 1.0,
                    'reactionsequence': u'PARENT',
                    'smiles': u'Oc1ccccc1O'
                },{
                    'isquery': True, 'level': 0, 'metid': 352, 'mol': u"Molfile of dihydroxyphenyl-valerolactone",
                    'molformula': u"C11H12O4",
                    'nhits': None,
                    'nr_scans': 1,
                    'origin': u"dihydroxyphenyl-valerolactone",
                    'probability': 1, 'reactionsequence': u"PARENT",
                    'smiles': u"O=C1OC(Cc2ccc(O)c(O)c2)CC1"
                }]
            }
        )

    def test_scanid(self):
        response = self.job.metabolites(scanid=641)
        self.assertIn('score', response['rows'][0])
        self.assertEqual(response['total'], 1)

    def test_filteredon_nrscanseq(self):
        response = self.job.metabolites(filters=[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}])
        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscansgt(self):
        response = self.job.metabolites(filters=[{"type":"numeric","comparison":"lt","value":2,"field":"nr_scans"}])
        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscanslt(self):
        response = self.job.metabolites(filters=[{"type":"numeric","comparison":"gt","value":0,"field":"nr_scans"}])
        self.assertEqual(response['total'], 2)

    def test_filteredon_isquery(self):
        response = self.job.metabolites(filters=[{"type":"boolean","value":True,"field":"isquery"}])
        self.assertEqual(response['total'], 2)

    def test_filteredon_molformula(self):
        response = self.job.metabolites(filters=[{"type":"string","value":"C6","field":"molformula"}])
        self.assertEqual(response['total'], 1)

    def test_filteredon_level(self):
        response = self.job.metabolites(filters=[{"type":"list","value":[0,1,2],"field":"level"}])
        self.assertEqual(response['total'], 2)

    def test_filteredon_score(self):
        response = self.job.metabolites(scanid=641, filters=[{"type":"numeric","comparison":"eq","value":200,"field":"score"}])
        self.assertEqual(response['total'], 1)

    def test_filteredon_score_without_scan(self):
        from mmm.job import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(filters=[{"type":"numeric","comparison":"eq","value":200,"field":"score"}])

    def test_sort_probmet(self):
        response = self.job.metabolites(sorts=[{"property":"probability","direction":"DESC"},{"property":"metid","direction":"ASC"}])
        self.assertEqual(response['total'], 2)

    def test_sort_nrscans(self):
        response = self.job.metabolites(sorts=[{"property":"nr_scans","direction":"DESC"}])
        self.assertEqual(response['total'], 2)

    def test_sort_score(self):
        response = self.job.metabolites(scanid=641, sorts=[{"property":"score","direction":"DESC"}])
        self.assertEqual(response['total'], 1)

    def test_sort_score_without_scan(self):
        from mmm.job import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(sorts=[{"property":"score","direction":"DESC"}])

class JobScansWithMetabolitesTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB())

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
        response = self.job.scansWithMetabolites(filters=[{"type":"string","value":"C6","field":"molformula"}])
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_nrscans(self):
        response = self.job.scansWithMetabolites(filters=[{"type":"numeric","value":"1", "comparison":"eq","field":"nr_scans"}])
        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_score(self):
        response = self.job.scansWithMetabolites(filters=[{"type":"numeric","value":"200", "comparison":"eq","field":"score"}])
        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317}
        ])

class jobMSpectraTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB())

    def test_scanonly(self):
        self.assertEqual(
            self.job.mspectra(641),
            {
                'peaks': [
                    {'intensity': 345608.65625, 'mz': 109.0295639038086},
                    {'intensity': 807576.625, 'mz': 305.033508300781}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': { 'id': None, 'mz': None }
            }
        )

    def test_withmslevel(self):
        self.assertEqual(
            self.job.mspectra(641, 1),
            {
                'peaks': [
                    {'intensity': 345608.65625, 'mz': 109.0295639038086},
                    {'intensity': 807576.625, 'mz': 305.033508300781}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': { 'id': None, 'mz': None }
            }
        )

    def test_notfound(self):
        from mmm.job import ScanNotFound
        with self.assertRaises(ScanNotFound):
            self.job.mspectra(123)

    def test_lvl2scan(self):
        response = self.job.mspectra(871, 2)
        self.assertEqual(response, {
            'peaks': [
                {'intensity': 211603.046875, 'mz': 123.04508972168},
                {'intensity': 279010.28125, 'mz': 163.076232910156}
            ],
            'cutoff': 139505.0,
            'mslevel': 2,
            'precursor': { 'id': 870, 'mz': 207.0663147 }
        })

class jobFragmentsTestCase(unittest.TestCase):
    def setUp(self):
        import uuid
        self.job = Job(uuid.uuid1(), initTestingDB())

    def test_metabolitewithoutfragments(self):
        response = self.job.fragments(metid=72, scanid=641)
        self.assertEqual(response, {
            'children': {
                'atoms': u'0,1,2,3,4,5,6,7',
                'children': [],
                'deltah': -1.0,
                'expanded': True,
                'fragid': 948,
                'leaf': True,
                'mass': 110.0367794368,
                'metid': 72,
                'mol': u'Molfile',
                'mslevel': 1,
                'mz': 109.0296783,
                'scanid': 641,
                'score': 200.0
            }, 'expanded': True
        })

    def test_metabolitewithfragments(self):
        response = self.job.fragments(metid=352, scanid=870)
        self.assertEqual(response, {
            'children': {
                'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14',
                'children': [{
                    'atoms': "6,7,8,9,10,11,12,13,14",
                    'deltah': 0,
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
                },{
                    'atoms': "3,4,5,6,7,8,9,10,11,12,13,14",
                    'deltah': -1,
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
                'expanded': True,
                'fragid': 1707,
                'leaf': False,
                'mass': 208.0735588736,
                'metid': 352,
                'mol': u'Molfile of dihydroxyphenyl-valerolactone',
                'mslevel': 1,
                'mz': 207.0663147,
                'scanid': 870,
                'score': 100
            }, 'expanded': True
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
            'score': 4
        }])

    def test_badfragment(self):
        from mmm.job import FragmentNotFound
        with self.assertRaises(FragmentNotFound):
            self.job.fragments(metid=70002, scanid=641)
