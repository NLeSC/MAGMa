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
            user = User(u'bob', u'Bob Example',
                        u'bob@example.com', 'mypassword')
            DBSession().add(user)
            self.job = jf.fromScratch(u'bob')
            self.jobid = self.job.id

    def do_login(self):
        params = {u'userid': u'bob', u'password': u'mypassword'}
        self.testapp.post('/login', params)

    def fake_jobid(self):
        """ Create job in self.root_dir filled with test db"""
        populateTestingDB(self.job.db.session)
        self.job.db.session.commit()
        return self.jobid

    def test_home(self):
        self.do_login()
        res = self.testapp.get('/', status=200)
        self.assertTrue('Submit' in res.body)

    def test_molecules(self):
        self.do_login()
        jobid = self.fake_jobid()

        res_url = '/results/' + str(jobid)
        res_url += '/molecules.json?limit=10&start=0'
        res = self.testapp.get(res_url, status=200)
        url1 = u'<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += u'?cid=289">CID: 289</a>'
        url2 = u'<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += u'?cid=152432">CID: 152432</a>'
        self.assertEqual(json.loads(res.body), {
            u'totalUnfiltered': 2,
            u'total': 2,
            u'scans': [{
                u'rt': 933.317,
                u'id': 641
            }, {
                u'rt': 1254.15,
                u'id': 870
            }],
            u'rows': [{
                u'molid': 72,
                u'predicted': False,
                u'mol': u'Molfile',
                u'formula': u'C6H6O2',
                u'nhits': 1,
                u'name': u'pyrocatechol',
                u'refscore': 1.0,
                u'reactionsequence': {
                                         u'reactantof': {
                                             u'esterase': {
                                                 u'nr': 2,
                                                 u'nrp': 1
                                             }
                                         }
                                     },
                u'smiles': u'C1=CC=C(C(=C1)O)O',
                u'inchikey14': u'YCIMNLLNPGFGHC',
                u'mim': 110.03677, u'logp': 1.231,
                u'assigned': False,
                u'reference': url1
            }, {
                u'predicted': False, 'molid': 352,
                u'mol': u"Molfile of dihydroxyphenyl-valerolactone",
                u'formula': u"C11H12O4",
                u'nhits': 1,
                u'name': u"dihydroxyphenyl-valerolactone",
                u'refscore': 1.0,
                u'reactionsequence': {
                                         u'productof': {
                                             u'theogallin': {
                                                 u'nr': 1,
                                                 u'nrp': 0
                                             }
                                         }
                                     },
                u'smiles': u"O=C1CCC(Cc2ccc(O)c(O)c2)O1",
                u'inchikey14': u'ZNXXWTPQHVLMQT',
                u'mim': 208.07355, u'logp': 2.763,
                u'assigned': False,
                u'reference': url2
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
        self.assertTrue('Submit' in res.body)
