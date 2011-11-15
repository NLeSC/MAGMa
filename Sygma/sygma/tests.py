import unittest
from pyramid.config import Configurator
from pyramid import testing

def _initTestingDB():
    from sqlalchemy import create_engine
    from sygma.models import initialize_sql, DBSession, Base
    engine = create_engine('sqlite://')
    initialize_sql(engine)
    Base.metadata.create_all(engine)
    session = DBSession
    _populateTestingDB(session)
    return session

def _populateTestingDB(session):
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
        basepeakmz=305.034, basepeakintensity=807577.0, totioncurrent=5957620,
    ))
    session.add_all([
        Peak(scanid=641, mz=305.033508300781, intensity=807576.625),
        Peak(scanid=641, mz=109.0295639038086, intensity=345608.65625)
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

class HomeView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def test_it(self):
        request = testing.DummyRequest()
        from sygma.views import home
        response = home(request)
        self.assertEqual(response, {})

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
        self.assertEqual(response['maxmslevel'], 1)

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
            'total': 1,
            'scans': [{
                'rt': 933.317,
                'id': 641
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
            }]
        })

    def test_filteredon_scanid(self):
        params = dict(start=0, limit=10, scanid=641)
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_filteredon_metid(self):
        params = dict(start=0, limit=10, metid=72)
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_filteredon_nrscanseq(self):
        params = dict(start=0, limit=10, filter='[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_filteredon_nrscansgt(self):
        params = dict(start=0, limit=10, filter='[{"type":"numeric","comparison":"lt","value":2,"field":"nr_scans"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_filteredon_nrscanslt(self):
        params = dict(start=0, limit=10, filter='[{"type":"numeric","comparison":"gt","value":0,"field":"nr_scans"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_filteredon_isquery(self):
        params = dict(start=0, limit=10, filter='[{"type":"boolean","value":true,"field":"isquery"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_filteredon_molformula(self):
        params = dict(start=0, limit=10, filter='[{"type":"string","value":"C6","field":"molformula"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_filteredon_level(self):
        params = dict(start=0, limit=10, filter='[{"type":"list","value":[0,1,2],"field":"level"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_sort_probmet(self):
        params = dict(start=0, limit=10, sort='[{"property":"probability","direction":"DESC"},{"property":"metid","direction":"ASC"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

    def test_sort_nrscans(self):
        params = dict(start=0, limit=10, sort='[{"property":"nr_scans","direction":"DESC"}]')
        response = self._callFUT(params)
        self.assertEqual(response['total'], 1)

class extracted_ion_chromatogram_QueryHelper(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, params):
        from sygma.views import extracted_ion_chromatogram
        return extracted_ion_chromatogram(params)

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
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

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
        self.assertEqual(self._callFUT(dict()), [{
            'id': 641,
            'rt': 933.317,
            'intensity': 807577.0
        }])

class ChromatogramHitsView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, params):
        from sygma.views import chromatogram_hits
        request = testing.DummyRequest(params=params)
        return chromatogram_hits(request)

    def test_allhits(self):
        params = dict()
        response = self._callFUT(params)
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_metid(self):
        params = dict(metid=72)
        response = self._callFUT(params)
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

class MSpectraView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_route('mspectra.json', '/mspectra/{id}.json')
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, id, params):
        from sygma.views import mspectrajson
        request = testing.DummyRequest(params=params)
        request.matchdict = {'id': id }
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


class MetaboliteScansView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_route('metabolite/scans.json','/metabolite/{id}/scans.json')
        self.session = _initTestingDB()

    def tearDown(self):
        self.session.remove()
        testing.tearDown()

    def _callFUT(self, id):
        from sygma.views import metabolitescans
        request = testing.DummyRequest()
        request.matchdict = {'id': id }
        return metabolitescans(request)

    def test_it(self):
        response = self._callFUT(72)
        self.assertEqual(response, {
            'chromatogram': [{
                'rt': 933.317,
                'intensity': 345608.65625
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

    def test_it(self):
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
            'total': 1,
            'scans': [{
                'rt': 933.317,
                'id': 641
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
            }]
        })