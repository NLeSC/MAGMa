#!/usr/bin/env python
"""Manage MAGMa web user accounts
"""

import argparse,sys
from transaction import commit
from paste.deploy import appconfig
from sqlalchemy import engine_from_config
from magmaweb.user import init_user_db, User
from magmaweb.job import JobFactory

class MagmaUserCommand(object):

    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('-c', '--config', help="magma web config file (default: production.ini)", default="production.ini", type=str)
        sp = self.parser.add_subparsers()

        sc_add = sp.add_parser("add", help=self.add.__doc__, description=self.add.__doc__)
        sc_add.add_argument('user', help="user id (default: %(default)s)", default=None,type=str)
        sc_add.add_argument('name', help="name (default: %(default)s)", default=None,type=str)
        sc_add.add_argument('email', help="e-mail address (default: %(default)s)", default=None,type=str)
        sc_add.add_argument('password', help="password (default: %(default)s)", default=None,type=str)
        sc_add.set_defaults(func=self.add)

        sc_update = sp.add_parser("update", help=self.update.__doc__, description=self.update.__doc__)
        sc_update.add_argument('-u', '--user', help="change user id (default: %(default)s)", default=None,type=str)
        sc_update.add_argument('-n', '--name', help="change name (default: %(default)s)", default=None,type=str)
        sc_update.add_argument('-e', '--email', help="change email (default: %(default)s)", default=None,type=str)
        sc_update.add_argument('-p', '--password', help="change password (default: %(default)s)", default=None,type=str)
        sc_update.add_argument('userid', help="user id (default: %(default)s)", default=None,type=str)
        sc_update.set_defaults(func=self.update)

        sc_remove = sp.add_parser("remove", help=self.remove.__doc__, description=self.remove.__doc__)
        sc_remove.add_argument('user', help="user id (default: %(default)s)", default=None,type=str)
        sc_remove.set_defaults(func=self.remove)

    def add(self, args):
        "Add new user"
        user = User(unicode(args.user), unicode(args.name), unicode(args.email), args.password)
        User.add(user)
        commit()

    def update(self, args):
        "Update user data"
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
        user=User.by_id(args.user)
        factory = JobFactory(args.datafolder+'jobs')
        for jobmeta in user.jobs:
            print jobmeta.jobid
            factory.fromId(jobmeta.jobid).delete()
        User.delete(user)
        commit()

    def run(self, argv):
        self.args = self.parser.parse_args(argv)

        self.config = appconfig('config:' + self.args.config, relative_to=)
        engine = engine_from_config(self.config)
        init_user_db(engine)

        print self.config
        exit(4)

        self.args.func(self.args)

def main(argv=sys.argv[1:]):
    command = MagmaUserCommand()
    return command.run(argv)
