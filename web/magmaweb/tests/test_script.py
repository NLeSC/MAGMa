from argparse import Namespace
from unittest import TestCase
from uuid import UUID
from mock import MagicMock
from transaction import commit
from magmaweb.job import JobFactory
from magmaweb.script import MagmaCommand
from magmaweb.user import DBSession, User, JobMeta
from .test_user import init_user_db, destroy_user_db


class FakeJob(object):
    id = '37dc6b15-2013-429c-98b7-f058bcf0c274'
    meta = JobMeta(UUID('37dc6b15-2013-429c-98b7-f058bcf0c274'), u'someone')


class TestMagmaSubCommands(TestCase):

    def setUp(self):
        self.command = MagmaCommand()
        init_user_db()
        self.command.job_factory = JobFactory('/dev/null')

    def tearDown(self):
        destroy_user_db()

    def testAddUser(self):
        args = Namespace(
            user=u'someone', name=u'somename', email=u'someemail', password='secret')

        self.command.add(args)

        self.assertEqual(DBSession().query(User.userid).all(), [(u'someone',)])

    def testUpdateUser_changeDisplayname(self):
        self.command.add(Namespace(
            user=u'someone', name=u'somename', email=u'someemail', password='secret'))

        args = Namespace(userid=u'someone', name=u'somename2')
        self.command.update(args)

        self.assertEqual(
            DBSession().query(User.displayname).all(), [(u'somename2',)])

    def testUpdateUser_changeEmail(self):
        self.command.add(Namespace(
            user=u'someone', name=u'somename', email=u'someemail', password='secret'))

        args = Namespace(userid=u'someone', email=u'someemail2')
        self.command.update(args)

        self.assertEqual(
            DBSession().query(User.email).all(), [(u'someemail2',)])

    def testRemoveUser(self):
        self.command.add(Namespace(
            user=u'someone', name=u'somename', email=u'someemail', password='secret'))

        args = Namespace(user=u'someone')
        self.command.remove(args)

        self.assertEqual(DBSession().query(User).all(), [])

    def testChangeJobOwner(self):
        user1 = User(u'someone', u'somename', u'someemail', 'secret')
        User.add(user1)
        user2 = User(u'someone2', u'somename2', u'someemail2', 'secret2')
        User.add(user2)
        job = FakeJob()
        self.command.job_factory.fromId = MagicMock(
            name='fromId', return_value=job)
        JobMeta.add(job.meta)
        commit()

        args = Namespace(
            job=u'37dc6b15-2013-429c-98b7-f058bcf0c274', user=u'someone2')
        self.command.owner(args)

        self.assertEqual(
            DBSession().query(JobMeta.owner).all(), [(u'someone2',)])

    def testImportJob(self):
        self.command.job_factory.fromDb = MagicMock(
            name='fromDb', return_value=FakeJob())

        args = Namespace(dbfile='results.db', owner=u'someone')
        self.command.importjob(args)

        self.command.job_factory.fromDb.assert_called_once_with(
            'results.db', u'someone')
