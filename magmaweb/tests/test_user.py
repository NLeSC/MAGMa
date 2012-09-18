import unittest
from mock import Mock
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
        j = user.Job('12345', 'My job')
        ju = user.JobUser('me', '12345', 'owner')
        self.session.add_all([u, j, ju])

    def tearDown(self):
        destroy_user_db()
        testing.tearDown()

    def test_getJob(self):
        from magmaweb.job import Job
        mjob = Mock(Job)
        jif = user.JobIdFactory(self.request)
        jif.job_factory.fromId = Mock(return_value=mjob)

        job = jif['12345']

        jif.job_factory.fromId.assert_called_once_with('12345')
        self.assertEqual(job, mjob)
        self.assertEqual(job.__parent__, jif)
        self.assertListEqual(job.__acl__, [(Allow, 'me', 'run')])

    def test_getJobNotFound(self):
        from magmaweb.job import JobNotFound
        jif = user.JobIdFactory(self.request)
        jif.job_factory.fromId = Mock(side_effect=JobNotFound('67890'))

        with self.assertRaises(HTTPNotFound):
            jif['67890']

    def make_authenticateduser_job_owner(self):
        jif = user.JobIdFactory(self.request)

        jif.make_authenticateduser_job_owner('2468')

        # look in db if owner set
        ju = self.session.query(user.JobUser).filter(user.JobUser.jobid=='2468').filter(user.JobUser.userid=='me').one()
        self.assertEqual(ju.role, 'owner')

class TestGroupFinder(unittest.TestCase):
    def setUp(self):
        init_user_db()
        self.session = user.DBSession()
        u = user.User('me', 'My', 'myself')
        ur = user.UserRole('me', 'admin')
        self.session.add_all([u, ur])

    def tearDown(self):
        destroy_user_db()

    def test_it(self):
        request = testing.DummyRequest()

        groups = user.groupfinder('me', request)

        self.assertListEqual(groups, ['admin'])
