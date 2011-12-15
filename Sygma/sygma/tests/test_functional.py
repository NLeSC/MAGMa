import unittest
from mock import patch

class FunctionalTests(unittest.TestCase):
    def setUp(self):
        from sygma import main
        settings = { 'mako.directories': 'sygma:templates'}
        app = main({}, **settings)
        from webtest import TestApp
        self.testapp = TestApp(app)

    def tearDown(self):
        del self.testapp

    def test_home(self):
        res = self.testapp.get('/', status=200)
        self.assertTrue('Homepage' in res.body)

    @patch('sygma.views.fetch_job')
    def test_metabolites(self, mocked_fetch_job):
        import uuid
        from test_job import initTestingDB
        from sygma.job import Job
        job = Job(uuid.uuid1(), initTestingDB())
        mocked_fetch_job.return_value = job

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