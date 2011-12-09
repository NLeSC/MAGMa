import unittest
from pyramid.config import Configurator
from pyramid import testing

def _initTestingDB(url = 'sqlite://'):
    """Creates testing db and populates with test data"""
    from sqlalchemy import create_engine
    from sygma.models import initialize_sql, DBSession, Base
    engine = create_engine(url) # default in memory db
    initialize_sql(engine)
    Base.metadata.create_all(engine)
    session = DBSession
    _populateTestingDB(session)
    return session

def _populateTestingDB(session):
    """Polulates test db with data

    Adds 1 metabolite with one fragment.
    Adds 1 metabolite with one fragment which has 2 child fragments of which one has another child fragment.
    Run, Scan and Peak are filled to create a working run.

    session
        session connection to db
    """
    from sygma.models import Metabolite, Scan, Peak, Fragment, Run
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

class HomeView(unittest.TestCase):
    def setUp(self):
        import tempfile
        settings = { 'jobrootdir': tempfile.mkdtemp() }
        self.config = testing.setUp(settings=settings)
        self.config.add_route('results', '/results')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.config.registry.settings['jobrootdir'])
        testing.tearDown()

    def test_get(self):
        request = testing.DummyRequest()
        from sygma.views import home
        response = home(request)
        self.assertEqual(response, {})

    def test_post(self):
        import os
        # post requires a sqlite db as file upload
        qdb = os.tmpnam()
        dbsession = _initTestingDB('sqlite:///'+qdb)
        dbsession.flush()
        dbsession.close()
        dbsession.remove()

        class FileUpload:
            pass

        post = { 'db': FileUpload() }
        post['db'].file = open(qdb,'r') # reopen tmpfile in readmode
        request = testing.DummyRequest(post=post)

        from sygma.views import home
        response = home(request)

        self.assertEqual(response.code, 302, 'Redirect')
        self.assertRegexpMatches(response.headers['Location'], '/results$', 'Redirected to results action')
        self.assertIn('id', request.session, 'Session has id')
        self.assertRegexpMatches(request.session['dbname'], 'results\.db$', 'Session has dbname')
        import filecmp
        rdbname = os.path.join(self.config.registry.settings['jobrootdir'], str(request.session['id']), 'results.db')
        self.assertTrue(filecmp.cmp(qdb, rdbname ), 'jobdir/<id>/results.db should be same as mocked input')
        os.remove(qdb)

class ResultsView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def test_it(self):
        request = testing.DummyRequest()
        from sygma.views import results
        response = results(request)
        self.assertEqual(response['run'].ms_intensity_cutoff, 200000.0)
        self.assertEqual(response['maxmslevel'], 3)

class MetabolitesView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, params):
        from sygma.views import metabolitesjson
        request = testing.DummyRequest(params=params)
        return metabolitesjson(request)

    def test_default(self):
        params = dict(start=0, limit=10)
        response = self._callFUT(params)
        self.assertEqual(response, {
            'total': 2,
            'scans': [{
                'rt': 933.317,
                'id': 641
            },{
               'rt': 1254.15,
               'id': 870
            }],
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
                'isquery': True, 'level': 0, 'metid': 352, 'mol': "Molfile of dihydroxyphenyl-valerolactone",
                'molformula': "C11H12O4",
                'nhits': None,
                'nr_scans': 1,
                'origin': "dihydroxyphenyl-valerolactone",
                'probability': 1, 'reactionsequence': "PARENT",
                'smiles': "O=C1OC(Cc2ccc(O)c(O)c2)CC1"
            }]
        })

    def test_filteredon_scanid(self):
        params = dict(start=0, limit=10, scanid=641)
        response = self._callFUT(params)
        self.assertIn('score', response['rows'][0])
        self.assertEqual(response['total'], 1)

    def test_filteredon_nrscanseq(self):
        params = dict(start=0, limit=10, filter='[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscansgt(self):
        params = dict(start=0, limit=10, filter='[{"type":"numeric","comparison":"lt","value":2,"field":"nr_scans"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscanslt(self):
        params = dict(start=0, limit=10, filter='[{"type":"numeric","comparison":"gt","value":0,"field":"nr_scans"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 2)

    def test_filteredon_isquery(self):
        params = dict(start=0, limit=10, filter='[{"type":"boolean","value":true,"field":"isquery"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 2)

    def test_filteredon_molformula(self):
        params = dict(start=0, limit=10, filter='[{"type":"string","value":"C6","field":"molformula"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_filteredon_level(self):
        params = dict(start=0, limit=10, filter='[{"type":"list","value":[0,1,2],"field":"level"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 2)

    def test_filteredon_score(self):
        params = dict(start=0, limit=10, scanid=641, filter='[{"type":"numeric","comparison":"eq","value":200,"field":"score"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_filteredon_score_without_scan(self):
        from sygma.views import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            params = dict(start=0, limit=10, filter='[{"type":"numeric","comparison":"eq","value":200,"field":"score"}]')
            response = self._callFUT(params)

    def test_sort_probmet(self):
        params = dict(start=0, limit=10, sort='[{"property":"probability","direction":"DESC"},{"property":"metid","direction":"ASC"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 2)

    def test_sort_nrscans(self):
        params = dict(start=0, limit=10, sort='[{"property":"nr_scans","direction":"DESC"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 2)

    def test_sort_score(self):
        params = dict(start=0, limit=10, scanid=641, sort='[{"property":"score","direction":"DESC"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_sort_score_without_scan(self):
        from sygma.views import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            params = dict(start=0, limit=10, sort='[{"property":"score","direction":"DESC"}]')
            response = self._callFUT(params)

class extracted_ion_chromatogram_QueryHelper(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, params):
        from sygma.views import filteredscans
        return filteredscans(testing.DummyRequest(), params)

    def test_metid(self):
        params = dict(metid=72)
        response = self._callFUT(params)
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_scanid(self):
        params = dict(scanid=641)
        response = self._callFUT(params)
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_molformula(self):
        params = dict(filter='[{"type":"string","value":"C6","field":"molformula"}]')
        response = self._callFUT(params)
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_nrscans(self):
        params = dict(filter='[{"type":"numeric","value":"1", "comparison":"eq","field":"nr_scans"}]')
        response = self._callFUT(params)
        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_score(self):
        params = dict(filter='[{"type":"numeric","value":"200", "comparison":"eq","field":"score"}]')
        response = self._callFUT(params)
        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317}
        ])

class ChromatogramView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, params):
        from sygma.views import chromatogramjson
        request = testing.DummyRequest(params=params)
        return chromatogramjson(request)

    def test_it(self):
        self.assertEqual(self._callFUT(dict()), [
            { 'id': 641, 'rt': 933.317, 'intensity': 807577.0 },
            { 'id': 870, 'rt': 1254.15, 'intensity': 1972180.0 }
        ])

class MSpectraView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_route('mspectra.json', '/mspectra/{scanid}.json')
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, id, params):
        from sygma.views import mspectrajson
        request = testing.DummyRequest(params=params)
        request.matchdict = {'scanid': id }
        return mspectrajson(request)

    def test_scanonly(self):
        response = self._callFUT(641,dict())
        self.assertEqual(response, {
            'peaks': [
                {'intensity': 345608.65625, 'mz': 109.0295639038086},
                {'intensity': 807576.625, 'mz': 305.033508300781}
            ],
            'cutoff': 200000.0,
            'mslevel': 1,
            'precursor': { 'id': None, 'mz': None }
        })

    def test_mslevel(self):
        response = self._callFUT(641,dict(mslevel=1))
        self.assertEqual(response, {
            'peaks': [
                {'intensity': 345608.65625, 'mz': 109.0295639038086},
                {'intensity': 807576.625, 'mz': 305.033508300781}
            ],
            'cutoff': 200000.0,
            'mslevel': 1,
            'precursor': { 'id': None, 'mz': None }
        })

    def test_notfound(self):
        self.assertEqual(self._callFUT(123, dict()).status, '404 Not Found')

    def test_lvl2scan(self):
        response = self._callFUT(871,dict())
        self.assertEqual(response, {
            'peaks': [
                {'intensity': 211603.046875, 'mz': 123.04508972168},
                {'intensity': 279010.28125, 'mz': 163.076232910156}
            ],
            'cutoff': 139505.0,
            'mslevel': 2,
            'precursor': { 'id': 870, 'mz': 207.0663147 }
        })


class MetaboliteScansView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_route('extractedionchromatogram.json','/extractedionchromatogram/{metid}.json')
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, id):
        from sygma.views import extractedionchromatogram
        request = testing.DummyRequest()
        request.matchdict = {'metid': id }
        return extractedionchromatogram(request)

    def test_it(self):
        response = self._callFUT(72)
        self.assertEqual(response, {
            'chromatogram': [{
                'rt': 933.317,
                'intensity': 345608.65625
            },{
                'rt': 1254.15,
                'intensity': 0
            }],
            'scans': [{
                'rt': 933.317,
                'id': 641
            }]
        })

    def test_badmetid(self):
        self.assertEqual(self._callFUT(123), {
            'chromatogram': [],
            'scans': []
        })

class FragmentsView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_route('fragments.json', '/fragments/{scanid}/{metid}.json')
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, metid, scanid, params):
        from sygma.views import fragments
        request = testing.DummyRequest(params=params)
        request.matchdict = {'metid': metid, 'scanid': scanid }
        return fragments(request)

    def test_metabolitewithoutfragments(self):
        response = self._callFUT(72, 641, dict(node=''))
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
        response = self._callFUT(352, 870, dict(node=''))
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
        response = self._callFUT(352, 870, dict(node=1709))
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

    def test_badmetabolite(self):
        response = self._callFUT(70002, 641, dict(node=''))
        self.assertEqual(response.status, '404 Not Found')

class FunctionalTests(unittest.TestCase):

    def setUp(self):
        from sygma import main
        settings = { 'sqlalchemy.url': 'sqlite:///:memory:', 'mako.directories': 'sygma:templates'}
        app = main({}, **settings)
        from sqlalchemy import engine_from_config
        from sygma.models import initialize_sql, DBSession, Base
        engine = engine_from_config(settings, 'sqlalchemy.')
        initialize_sql(engine)
        Base.metadata.create_all(engine)
        _populateTestingDB(DBSession)
        from webtest import TestApp
        self.testapp = TestApp(app)

    def tearDown(self):
        del self.testapp
        from sygma.models import DBSession
        DBSession.remove()

    def test_home(self):
        res = self.testapp.get('/', status=200)
        self.assertTrue('Homepage' in res.body)

    def test_metabolites(self):
        res = self.testapp.get('/metabolites.json?limit=10&start=0', status=200)
        import simplejson as json
        self.assertEqual(json.loads(res.body),{
            'total': 2,
            'scans': [{
                'rt': 933.317,
                'id': 641
            },{
               'rt': 1254.15,
               'id': 870
            }],
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
                'isquery': True, 'level': 0, 'metid': 352, 'mol': "Molfile of dihydroxyphenyl-valerolactone",
                'molformula': "C11H12O4",
                'nhits': None,
                'nr_scans': 1,
                'origin': "dihydroxyphenyl-valerolactone",
                'probability': 1, 'reactionsequence': "PARENT",
                'smiles': "O=C1OC(Cc2ccc(O)c(O)c2)CC1"
            }]
        })