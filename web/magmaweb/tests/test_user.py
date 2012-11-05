import unittest
import uuid
from mock import Mock, patch
from sqlalchemy import create_engine
from pyramid import testing
from pyramid.httpexceptions import HTTPNotFound
from pyramid.security import Allow, Deny, Everyone, ALL_PERMISSIONS, Authenticated
from magmaweb import user

def init_user_db():
    engine = create_engine('sqlite:///:memory:')
    user.DBSession.configure(bind=engine)
    user.Base.metadata.create_all(engine) #@UndefinedVariable

def destroy_user_db():
    user.DBSession.remove()

class TestRootFactory(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_static_view('static', 'magmaweb:static')
        self.request = testing.DummyRequest()
        self.request.registry.settings = {'extjsroot': 'extjsroot'}

    def tearDown(self):
        destroy_user_db()
        testing.tearDown()

    def test_acl(self):
        rf = user.RootFactory(self.request)

        expected_acl = [
            (Allow, Authenticated, 'view'),
            (Deny, Everyone, ALL_PERMISSIONS)
        ]
        self.assertEqual(expected_acl, rf.__acl__)

    def test_extjsroot(self):
        rf = user.RootFactory(self.request)

        self.assertEqual(rf.request.extjsroot, 'http://example.com/static/extjsroot')

    def test_user_nologin(self):
        rf = user.RootFactory(self.request)

        self.assertIsNone(rf.request.user)

    @patch('magmaweb.user.unauthenticated_userid')
    def test_user_notindb(self, uau):
        uau.return_value = 'bob'
        init_user_db()

        rf = user.RootFactory(self.request)

        self.assertIsNone(rf.request.user)

    @patch('magmaweb.user.unauthenticated_userid')
    def test_user(self, uau):
        uau.return_value = 'bob'
        init_user_db()
        u = user.User('bob', 'Bobs name', 'bob@example.com')
        user.DBSession().add(u)

        rf = user.RootFactory(self.request)

        self.assertEqual(rf.request.user, u)


class TestJobIdFactory(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_static_view('static', 'magmaweb:static')
        self.request = testing.DummyRequest()
        self.request.registry.settings = {
                                          'extjsroot': 'extjsroot',
                                          'jobfactory.root_dir': '/somedir',
        }
        init_user_db()
        self.session = user.DBSession()
        u = user.User('me', 'My', 'myself')
        self.session.add(u)
        job_id = uuid.UUID('11111111-1111-1111-1111-111111111111')
        j = user.JobMeta(job_id, 'My job', owner='me')
        self.session.add(j)

    def tearDown(self):
        destroy_user_db()
        testing.tearDown()

    def test_getJob(self):
        from magmaweb.job import Job
        mjob = Mock(Job)
        mjob.owner = 'bob'
        jif = user.JobIdFactory(self.request)
        jif.job_factory.fromId = Mock(return_value=mjob)
        job_id = uuid.UUID('11111111-1111-1111-1111-111111111111')

        job = jif[str(job_id)]

        jif.job_factory.fromId.assert_called_once_with(job_id)
        self.assertEqual(job, mjob)
        self.assertEqual(job.__parent__, jif)
        self.assertEqual(job.__acl__, [(Allow, 'bob', 'run')])

    def test_getJobNotFound(self):
        from magmaweb.job import JobNotFound
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        jif = user.JobIdFactory(self.request)
        notfound = JobNotFound('Job not found', jobid)
        jif.job_factory.fromId = Mock(side_effect=notfound)

        with self.assertRaises(HTTPNotFound):
            jif[str(jobid)]

class TestUtils(unittest.TestCase):
    def setUp(self):
        init_user_db()
        self.session = user.DBSession()
        u = user.User('me', 'My', 'myself')
        self.session.add(u)
        job_id = uuid.UUID('11111111-1111-1111-1111-111111111111')
        j = user.JobMeta(job_id, 'My job', owner='me')
        self.session.add(j)

    def tearDown(self):
        destroy_user_db()

    def test_jobs_of_user(self):
        u = self.session.query(user.User).get('me')
        jobs = u.jobs
        self.assertEqual(jobs[0].description, 'My job')
