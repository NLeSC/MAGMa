#!/usr/bin/env python
"""Manage MAGMa web user accounts
"""

import argparse,sys
from transaction import commit
from sqlalchemy import create_engine
from magmaweb.user import init_user_db, User
from magmaweb.job import JobFactory

class MagmaUserCommand(object):

    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('config', help="magma web config file (default: production.ini)", default="production.ini", type=str)

        sc_add = self.parser.add_subparsers("add", help=self.add.__doc__, description=self.add.__doc__)
        sc_add.add_argument('user', help="user id (default: %(default)s)", default=None,type=str)
        sc_add.add_argument('name', help="name (default: %(default)s)", default=None,type=str)
        sc_add.add_argument('email', help="e-mail address (default: %(default)s)", default=None,type=str)
        sc_add.add_argument('password', help="password (default: %(default)s)", default=None,type=str)
        sc_add.add_argument('datafolder', help="folder containing users.db and jobs folder (default: %(default)s)", default='./',type=str)
        sc_add.set_defaults(func=self.add)

        sc_update = self.parser.add_subparsers("update", help=self.update.__doc__, description=self.update.__doc__)
        sc_update.add_argument('-u', '--user', help="change user id (default: %(default)s)", default=None,type=str)
        sc_update.add_argument('-n', '--name', help="change name (default: %(default)s)", default=None,type=str)
        sc_update.add_argument('-e', '--email', help="change email (default: %(default)s)", default=None,type=str)
        sc_update.add_argument('-p', '--password', help="change password (default: %(default)s)", default=None,type=str)
        sc_update.add_argument('userid', help="user id (default: %(default)s)", default=None,type=str)
        sc_update.add_argument('datafolder', help="folder containing users.db and jobs folder (default: %(default)s)", default='./',type=str)
        sc_update.set_defaults(func=self.update)

        sc_remove = self.parser.add_subparsers("remove", help=self.remove.__doc__, description=self.remove.__doc__)
        sc_remove.add_argument('user', help="user id (default: %(default)s)", default=None,type=str)
        sc_remove.add_argument('datafolder', help="folder containing users.db and jobs folder (default: %(default)s)", default='./',type=str)
        sc_remove.set_defaults(func=self.remove)

    def add(self, args):
        "Add new user"
        s='sqlite:///'+args.datafolder+'users.db'
        engine = create_engine(s)
        init_user_db(engine)

        user = User(unicode(args.user),
        unicode(args.name),
        unicode(args.email),
        args.password)
        User.add(user)
        commit()

    def update(self, args):
        "Update user data"
        s='sqlite:///'+args.datafolder+'users.db'
        engine = create_engine(s)
        init_user_db(engine)
        user=User.by_id(args.userid)
        if args.user != None:
            factory = JobFactory(args.datafolder+'jobs')
        for jobmeta in user.jobs:
            print jobmeta.jobid
        factory.fromId(jobmeta.jobid).set_owner(args.user)
        user.userid = unicode(args.user)
        if args.name != None:
            user.displayname = unicode(args.name)
        if args.email != None:
            user.email = unicode(args.email)
        if args.password != None:
            user.password = args.password
        User.add(user)
        commit()

    def remove(self, args):
        "Remove user and his/her jobs"
        s='sqlite:///'+args.datafolder+'users.db'
        engine = create_engine(s)
        init_user_db(engine)
        user=User.by_id(args.user)
        factory = JobFactory(args.datafolder+'jobs')
        for jobmeta in user.jobs:
            print jobmeta.jobid
        factory.fromId(jobmeta.jobid).delete()
        User.delete(user)
        commit()

    def run(self, argv):
        self.args = self.parser.parse_args(argv)
        self.args.func(self.args)

def main(argv=sys.argv):
    command = MagmaUserCommand()
    return command.run(argv)
