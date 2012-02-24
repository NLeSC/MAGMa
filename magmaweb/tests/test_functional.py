import unittest
from mock import patch

class FunctionalTests(unittest.TestCase):
    def setUp(self):
        from magmaweb import main
        settings = { 'mako.directories': 'magmaweb:templates', 'extjsroot': 'ext' }
        app = main({}, **settings)
        from webtest import TestApp
        self.testapp = TestApp(app)

    def tearDown(self):
        del self.testapp

    def test_home(self):
        res = self.testapp.get('/', status=200)
        self.assertTrue('Homepage' in res.body)

    @patch('magmaweb.views.fetch_job')
    def test_metabolites(self, mocked_fetch_job):
        import uuid, tempfile, shutil
        from test_job import initTestingDB
        from magmaweb.job import Job
        jobdir = tempfile.mkdtemp()
        job = Job(uuid.uuid1(), initTestingDB(), jobdir)
        mocked_fetch_job.return_value = job

        res = self.testapp.get('/results/'+str(job.id)+'/metabolites.json?limit=10&start=0', status=200)
        import json
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
                'smiles': u'Oc1ccccc1O',
                'mim': 110.03677, 'logp':1.231
            },{
                'isquery': True, 'level': 0, 'metid': 352, 'mol': "Molfile of dihydroxyphenyl-valerolactone",
                'molformula': "C11H12O4",
                'nhits': None,
                'nr_scans': 1,
                'origin': "dihydroxyphenyl-valerolactone",
                'probability': 1, 'reactionsequence': "PARENT",
                'smiles': "O=C1OC(Cc2ccc(O)c(O)c2)CC1",
                'mim': 208.07355, 'logp':2.763
            }]
        })

        shutil.rmtree(jobdir)