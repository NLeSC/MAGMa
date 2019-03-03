import tempfile
import unittest
import transaction
import json
from nose.plugins.attrib import attr
from webtest import TestApp
from magmaweb import main
from magmaweb.user import DBSession, User
from magmaweb.job import make_job_factory
from magmaweb.tests.test_job import populateTestingDB


@attr('functional')
class FunctionalTests(unittest.TestCase):
    settings = {}

    def setUp(self):
        self.root_dir = tempfile.mkdtemp()
        # default settings
        settings = {'jobfactory.root_dir': self.root_dir,
                    'mako.directories': 'magmaweb:templates',
                    'extjsroot': 'ext',
                    'sqlalchemy.url': 'sqlite:///:memory:',
                    'cookie.secret': 'aepeeV6aizaiph5Ae0Reimeequuluwoh',
                    'cookie.path': '/',
                    'monitor_user': 'jobmanager',
                    }
        settings.update(self.settings)
        self.settings = settings
        app = main({}, **self.settings)
        self.testapp = TestApp(app)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.root_dir)
        del self.testapp
        DBSession.remove()


class FunctionalPrivateTests(FunctionalTests):
    def setUp(self):
        FunctionalTests.setUp(self)
        # Setup owner of job
        jf = make_job_factory(self.settings)
        with transaction.manager:
            user = User('bob', 'Bob Example',
                        'bob@example.com', 'mypassword')
            DBSession().add(user)
            self.job = jf.fromScratch('bob')
            self.jobid = self.job.id

    def do_login(self):
        params = {'userid': 'bob', 'password': 'mypassword'}
        self.testapp.post('/login', params)

    def fake_jobid(self):
        """ Create job in self.root_dir filled with test db"""
        with transaction.manager:
            populateTestingDB(self.job.db.session)
        return self.jobid

    def test_home(self):
        self.do_login()
        res = self.testapp.get('/', status=200)
        self.assertTrue(b'Submit' in res.body)

    def test_molecules(self):
        self.do_login()
        jobid = self.fake_jobid()

        res_url = '/results/' + str(jobid)
        res_url += '/molecules.json?limit=10&start=0'
        res = self.testapp.get(res_url, status=200)
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += '?cid=152432">CID: 152432</a>'
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
                'molid': 72,
                'predicted': False,
                'mol': 'Molfile',
                'formula': 'C6H6O2',
                'nhits': 1,
                'name': 'pyrocatechol',
                'refscore': 1.0,
                'reactionsequence': {
                                         'reactantof': {
                                             'esterase': {
                                                 'nr': 2,
                                                 'nrp': 1
                                             }
                                         }
                                     },
                'smiles': 'C1=CC=C(C(=C1)O)O',
                'inchikey14': 'YCIMNLLNPGFGHC',
                'mim': 110.03677, 'logp': 1.231,
                'assigned': False,
                'reference': url1
            }, {
                'predicted': False, 'molid': 352,
                'mol': "Molfile of dihydroxyphenyl-valerolactone",
                'formula': "C11H12O4",
                'nhits': 1,
                'name': "dihydroxyphenyl-valerolactone",
                'refscore': 1.0,
                'reactionsequence': {
                                         'productof': {
                                             'theogallin': {
                                                 'nr': 1,
                                                 'nrp': 0
                                             }
                                         }
                                     },
                'smiles': "O=C1CCC(Cc2ccc(O)c(O)c2)O1",
                'inchikey14': 'ZNXXWTPQHVLMQT',
                'mim': 208.07355, 'logp': 2.763,
                'assigned': False,
                'reference': url2
            }]
        })

    def test_double_update_job(self):
        """Double update should not raise
        OperationalError: database is locked"""
        self.do_login()
        jobid = self.fake_jobid()
        req_url = '/results/' + str(jobid)
        req_body = json.dumps({"id": "bar",
                               "description": "New description",
                               "ms_filename": "F6789.mzxml",
                               "created_at": "1999-12-17T13:45:04",
                               "is_public": False,
                               })
        self.testapp.put(req_url, req_body)

        req_body2 = json.dumps({"id": "bar",
                                "description": "New description 2",
                                "ms_filename": "F6789.mzxml 2",
                                "created_at": "1999-12-17T13:45:04",
                                "is_public": False,
                                })

        self.testapp.put(req_url, req_body2)


class FunctionalPublicTests(FunctionalTests):
    settings = {'auto_register': True}

    def test_home(self):
        # Visiting / redirects to /login
        # which automatically registers/logins and redirects back to /
        res = self.testapp.get('/', status=302).follow(status=200)
        self.assertTrue(b'Submit' in res.body)
