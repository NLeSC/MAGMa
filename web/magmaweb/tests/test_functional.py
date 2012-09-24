import unittest

class FunctionalTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        from magmaweb import main
        self.root_dir = tempfile.mkdtemp()
        self.settings = {
                    'jobfactory.root_dir': self.root_dir,
                    'mako.directories': 'magmaweb:templates', 'extjsroot': 'ext'
                    }
        app = main({}, **self.settings)
        from webtest import TestApp
        self.testapp = TestApp(app)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.root_dir)
        del self.testapp

    def test_home(self):
        res = self.testapp.get('/', status=200)
        self.assertTrue('Homepage' in res.body)

    def fake_jobid(self):
        """ Create job in self.root_dir filled with test db"""
        from magmaweb.job import make_job_factory
        jf = make_job_factory(self.settings)
        job = jf.fromScratch()
        from test_job import populateTestingDB
        populateTestingDB(job.session)
        job.session.commit()
        return job.id

    def test_metabolites(self):
        jobid = self.fake_jobid()

        res = self.testapp.get('/results/'+str(jobid)+'/metabolites.json?limit=10&start=0', status=200)
        import json
        self.assertEqual(json.loads(res.body),{
            'totalUnfiltered': 2,
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
                'nhits': 1,
                'origin': u'pyrocatechol',
                'probability': 1.0,
                'reactionsequence': u'PARENT',
                'smiles': u'Oc1ccccc1O',
                'mim': 110.03677, 'logp':1.231,
                'assigned': False,
                'reference': '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=289">CID: 289</a>'
            },{
                'isquery': True, 'level': 0, 'metid': 352, 'mol': "Molfile of dihydroxyphenyl-valerolactone",
                'molformula': "C11H12O4",
                'nhits': None,
                'nhits': 1,
                'origin': "dihydroxyphenyl-valerolactone",
                'probability': 1, 'reactionsequence': "PARENT",
                'smiles': "O=C1OC(Cc2ccc(O)c(O)c2)CC1",
                'mim': 208.07355, 'logp':2.763,
                'assigned': False,
                'reference': '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=152432">CID: 152432</a>',
            }]
        })