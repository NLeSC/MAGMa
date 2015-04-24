#!/usr/bin/env python
"""Manage MAGMa web user accounts
"""

import argparse
import os
import sys
from transaction import commit
from paste.deploy import appconfig
from sqlalchemy import engine_from_config
from magmaweb.user import init_user_db, User, JobMeta
from magmaweb.job import make_job_factory


class MagmaCommand(object):

    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument(
            '-c', '--config',
            help="magma web config file (default: production.ini)",
            default="production.ini", type=str)
        sp = self.parser.add_subparsers()

        sc_add = sp.add_parser(
            "add", help=self.add.__doc__, description=self.add.__doc__)
        sc_add.add_argument(
            'user', help="user id (default: %(default)s)",
            default=None, type=unicode)
        sc_add.add_argument(
            'name', help="name (default: %(default)s)",
            default=None, type=unicode)
        sc_add.add_argument(
            'email', help="e-mail address (default: %(default)s)",
            default=None, type=unicode)
        sc_add.add_argument(
            'password', help="password (default: %(default)s)",
            default=None, type=str)
        sc_add.set_defaults(func=self.add)

        sc_update = sp.add_parser(
            "update", help=self.update.__doc__,
            description=self.update.__doc__)
        sc_update.add_argument(
            '-u', '--user', help="change user id (default: %(default)s)",
            default=argparse.SUPPRESS, type=unicode)
        sc_update.add_argument(
            '-n', '--name', help="change display name (default: %(default)s)",
            default=argparse.SUPPRESS, type=unicode)
        sc_update.add_argument(
            '-e', '--email', help="change email (default: %(default)s)",
            default=argparse.SUPPRESS, type=unicode)
        sc_update.add_argument(
            '-p', '--password', help="change password (default: %(default)s)",
            default=argparse.SUPPRESS, type=str)
        sc_update.add_argument(
            'userid', help="user id (default: %(default)s)",
            default=None, type=unicode)
        sc_update.set_defaults(func=self.update)

        sc_remove = sp.add_parser(
            "remove", help=self.remove.__doc__,
            description=self.remove.__doc__)
        sc_remove.add_argument(
            'user', help="user id (default: %(default)s)",
            default=None, type=unicode)
        sc_remove.set_defaults(func=self.remove)

        sc_owner = sp.add_parser(
            "owner", help=self.owner.__doc__, description=self.owner.__doc__)
        sc_owner.add_argument(
            'job', help="job identifier", default=None, type=unicode)
        sc_owner.add_argument(
            'user', help="user id", default=None, type=unicode)
        sc_owner.set_defaults(func=self.owner)

        sc_import = sp.add_parser(
            "importjob", help=self.importjob.__doc__,
            description=self.importjob.__doc__)
        sc_import.add_argument(
            'dbfile', help="job sqlite result db file",
            default=None, type=argparse.FileType('r'))
        sc_import.add_argument(
            'owner', help="user id", default=None, type=unicode)
        sc_import.set_defaults(func=self.importjob)

    def add(self, args):
        "Add new user"
        user = User(args.user, args.name, args.email, args.password)
        User.add(user)
        commit()

    def update(self, args):
        "Update user data"
        user = User.by_id(args.userid)
        if 'user' in args:
            user.userid = args.user
            for job in user.jobs:
                job.owner = args.user
                JobMeta.add(job)
        if 'name' in args:
            user.displayname = args.name
        if 'email' in args:
            user.email = args.email
        if 'password' in args:
            user.password = args.password
        User.add(user)
        commit()

    def remove(self, args):
        "Remove user and his/her jobs"
        user = User.by_id(args.user)
        for jobmeta in user.jobs:
            print jobmeta.jobid
            self.job_factory.fromId(jobmeta.jobid).delete()
        User.delete(user)
        commit()

    def owner(self, args):
        """Alter owner of job"""
        job = self.job_factory.fromId(args.job)
        job.meta.owner = args.user
        JobMeta.add(job.meta)
        commit()

    def importjob(self, args):
        """Import job results db"""
        dbfile = args.dbfile
        owner = args.owner
        job = self.job_factory.fromDb(dbfile, owner)
        print job.id
        commit()

    def configure(self, config_file):
        config_url = 'config:' + config_file
        cwd = os.getcwd()
        self.config = appconfig(config_url, 'MAGMaWeb', relative_to=cwd)
        engine = engine_from_config(self.config)
        init_user_db(engine)
        self.job_factory = make_job_factory(self.config)

    def run(self, argv):
        args = self.parser.parse_args(argv)

        self.configure(args.config)

        args.func(args)


def main(argv=sys.argv[1:]):
    command = MagmaCommand()
    return command.run(argv)
