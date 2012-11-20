import unittest
import uuid
import datetime
from mock import Mock, patch
from sqlalchemy import create_engine
from pyramid import testing
from pyramid.httpexceptions import HTTPNotFound
from pyramid.security import Allow, Deny, Everyone
from pyramid.security import ALL_PERMISSIONS, Authenticated
from magmaweb import user


def init_user_db():
    engine = create_engine('sqlite:///:memory:')
    user.init_user_db(engine, True, False)


def destroy_user_db():
    user.DBSession.remove()


class TestUser(unittest.TestCase):
    def test_construct(self):
        u = user.User('bob', 'Bob Smith', 'bob@smith.org')
        self.assertEqual(u.userid, 'bob')
        self.assertEqual(u.displayname, 'Bob Smith')
        self.assertEqual(u.email, 'bob@smith.org')

    @patch('bcrypt.hashpw')
    @patch('bcrypt.gensalt')
    def test_construct_with_password(self, salt, hashpw):
        salt.return_value = 'salty'

        u = user.User('bob', 'Bob Smith', 'bob@smith.org', 'mypassword')

        self.assertEqual(u.userid, 'bob')
        self.assertEqual(u.displayname, 'Bob Smith')
        self.assertEqual(u.email, 'bob@smith.org')
        salt.assert_called_with(12)
        hashpw.assert_called_with('mypassword', 'salty')

    def test_validate_password_correct(self):
        u = user.User('bob', 'Bob Smith', 'bob@smith.org', 'mypassword')
        self.assertTrue(u.validate_password('mypassword'))

    def test_validate_password_incorrect(self):
        u = user.User('bob', 'Bob Smith', 'bob@smith.org', 'mypassword')
        self.assertFalse(u.validate_password('otherpassword'))

    def test_repr(self):
        u = user.User('bob', 'Bob Smith', 'bob@smith.org', 'mypassword')
        self.assertEqual(repr(u),
                         "<User('bob', 'Bob Smith', 'bob@smith.org')>")

    def test_by_id(self):
        init_user_db()
        self.session = user.DBSession()
        u = user.User('bob', 'Bob Smith', 'bob@smith.org')
        self.session.add(u)

        u2 = user.User.by_id('bob')
        self.assertEqual(u, u2)

        destroy_user_db()

    def test_add(self):
        init_user_db()

        u = user.User('bob', 'Bob Smith', 'bob@smith.org')
        user.User.add(u)

        session = user.DBSession()
        u2 = session.query(user.User).get('bob')
        self.assertEqual(u, u2)

        destroy_user_db()

    def test_jobs(self):
        init_user_db()
        session = user.DBSession()
        u = user.User('bob', 'Bob Smith', 'bob@smith.org')
        session.add(u)
        job_id = uuid.UUID('11111111-1111-1111-1111-111111111111')
        j = user.JobMeta(job_id, 'bob')
        session.add(j)

        u2 = user.User.by_id('bob')  # force commit
        self.assertEqual(u2.jobs, [j])

        destroy_user_db()


class TestJobMeta(unittest.TestCase):
    @patch('datetime.datetime')
    def test_contruct_minimal(self, mock_dt):
        jid = uuid.UUID('986917b1-66a8-42c2-8f77-00be28793e58')
        created_at = datetime.datetime(2012, 11, 14, 10, 48, 26, 504478)
        mock_dt.utcnow.return_value = created_at

        j = user.JobMeta(jid, 'bob')

        self.assertEqual(j.jobid, jid)
        self.assertEqual(j.owner, 'bob')
        self.assertEqual(j.description, '')
        self.assertEqual(j.ms_filename, '')
        self.assertIsNone(j.parentjobid)
        self.assertEqual(j.state, 'STOPPED')
        self.assertEqual(j.created_at, created_at)

    def test_contstruct(self):
        jid = uuid.UUID('986917b1-66a8-42c2-8f77-00be28793e58')
        pid = uuid.UUID('83198655-b287-427f-af0d-c6bc1ca566d8')
        created_at = datetime.datetime(2012, 11, 14, 10, 48, 26, 504478)

        j = user.JobMeta(jid, 'bob', description='My desc',
                         parentjobid=pid, state='RUNNING',
                         ms_filename='F00346.mzxml',
                         created_at=created_at)

        self.assertEqual(j.jobid, jid)
        self.assertEqual(j.owner, 'bob')
        self.assertEqual(j.description, 'My desc')
        self.assertEqual(j.ms_filename, 'F00346.mzxml')
        self.assertEqual(j.parentjobid, pid)
        self.assertEqual(j.state, 'RUNNING')
        self.assertEqual(j.created_at, created_at)


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

        self.assertEqual(rf.request.extjsroot,
                         'http://example.com/static/extjsroot')

    def test_user_nologin(self):
        rf = user.RootFactory(self.request)

        self.assertIsNone(rf.request.user)

    @patch('magmaweb.user.authenticated_userid')
    def test_user_notindb(self, uau):
        uau.return_value = 'bob'
        init_user_db()

        rf = user.RootFactory(self.request)

        self.assertIsNone(rf.request.user)

    @patch('magmaweb.user.authenticated_userid')
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
        self.request.registry.settings = {'extjsroot': 'extjsroot',
                                          'jobfactory.root_dir': '/somedir',
                                          }
        init_user_db()
        self.session = user.DBSession()
        u = user.User('me', 'My', 'myself')
        self.session.add(u)
        job_id = uuid.UUID('11111111-1111-1111-1111-111111111111')
        j = user.JobMeta(job_id, 'me', description='My job')
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
