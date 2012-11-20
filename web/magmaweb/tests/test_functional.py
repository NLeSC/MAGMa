import tempfile
import unittest
import transaction
import json
from webtest import TestApp
from magmaweb import main
from magmaweb.user import DBSession, User
from magmaweb.job import make_job_factory
from test_job import populateTestingDB


class FunctionalTests(unittest.TestCase):
    def setUp(self):
        self.root_dir = tempfile.mkdtemp()
        self.settings = {
                         'jobfactory.root_dir': self.root_dir,
                         'mako.directories': 'magmaweb:templates',
                         'extjsroot': 'ext',
                         'sqlalchemy.url': 'sqlite:///:memory:',
                         'cookie.secret': 'aepeeV6aizaiph5Ae0Reimeequuluwoh',
                         'cookie.path': '/magma'
                         }
        # add REMOTE_ADDR to prevent ip auth
        # treating requests as coming from localhost
        env = dict(REMOTE_ADDR='1.2.3.4')
        app = main({}, **self.settings)
        self.testapp = TestApp(app, extra_environ=env)

        # Setup owner of job
        jf = make_job_factory(self.settings)
        with transaction.manager:
            user = User('bob', 'Bob Example', 'bob@example.com', 'mypassword')
            DBSession().add(user)
            self.job = jf.fromScratch('bob')
            self.jobid = self.job.id

    def tearDown(self):
        import shutil
        shutil.rmtree(self.root_dir)
        del self.testapp
        DBSession.remove()

    def test_home(self):
        res = self.testapp.get('/', status=200)
        self.assertTrue('Welcome' in res.body)

    def fake_jobid(self):
        """ Create job in self.root_dir filled with test db"""
        populateTestingDB(self.job.db.session)
        self.job.db.session.commit()
        return self.jobid

    def do_login(self):
        params = {'userid': 'bob', 'password': 'mypassword'}
        self.testapp.post('/login', params)

    def test_metabolites(self):
        self.do_login()
        jobid = self.fake_jobid()

        res_url = '/results/' + str(jobid)
        res_url += '/metabolites.json?limit=10&start=0'
        res = self.testapp.get(res_url, status=200)
        self.assertEqual(json.loads(res.body), {
            'totalUnfiltered': 2,
            'total': 2,
            'scans': [{
                'rt': 933.317,
                'id': 641
            }, {
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
                'reactionsequence': [u'PARENT'],
                'smiles': u'Oc1ccccc1O',
                'mim': 110.03677, 'logp':1.231,
                'assigned': False,
                'reference': '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=289">CID: 289</a>'
            }, {
                'isquery': True, 'level': 0, 'metid': 352,
                'mol': "Molfile of dihydroxyphenyl-valerolactone",
                'molformula': "C11H12O4",
                'nhits': None,
                'nhits': 1,
                'origin': "dihydroxyphenyl-valerolactone",
                'probability': 1,
                'reactionsequence': ["PARENT", "CHILD"],
                'smiles': "O=C1OC(Cc2ccc(O)c(O)c2)CC1",
                'mim': 208.07355, 'logp':2.763,
                'assigned': False,
                'reference': '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=152432">CID: 152432</a>',
            }]
        })
