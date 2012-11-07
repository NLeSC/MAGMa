import uuid
from sqlalchemy import Column
from sqlalchemy import Unicode
from sqlalchemy import ForeignKey
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.orm import (
    scoped_session,
    sessionmaker,
    )
import sqlalchemy.types as types
from zope.sqlalchemy import ZopeTransactionExtension #@UnresolvedImport
from pyramid.httpexceptions import HTTPNotFound
from pyramid.security import Allow, Deny, Everyone, ALL_PERMISSIONS
from pyramid.security import Authenticated, unauthenticated_userid
from magmaweb.job import make_job_factory, JobNotFound

DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))

Base = declarative_base()


class UUIDType(types.TypeDecorator):
    """ Stores UUID as Unicode string in database"""
    impl = Unicode

    def process_bind_param(self, value, dialect):
        if value is None:
            return None
        return str(value)

    def process_result_value(self, value, dialect):
        if value is None:
            return None
        return uuid.UUID(value)


class User(Base):
    """User list"""
    __tablename__= 'user'
    userid = Column(Unicode, primary_key=True)
    displayname = Column(Unicode)
    email = Column(Unicode)
    jobs = relationship("JobMeta", order_by="JobMeta.jobid")

    def __init__(self, userid, displayname, email):
        self.userid = userid
        self.displayname = displayname
        self.email = email


class JobMeta(Base):
    """Job metadata"""
    __tablename__ = 'job'
    jobid = Column(UUIDType, primary_key=True)
    description = Column(Unicode)  # Kept in sync with Run.description
    parentjobid = Column(UUIDType, ForeignKey('job.jobid'))
    state = Column(Unicode)  # queued/running/ready etc.
    owner = Column(Unicode, ForeignKey('user.userid'))
    children = relationship('JobMeta',
                            backref=backref('parent', remote_side=[jobid]))

    def __init__(self, jobid,
                 description='', parentjobid=None,
                 state=None, owner=None):
        self.jobid = jobid
        self.description = description
        self.parentjobid = parentjobid
        self.state = state
        self.owner = owner


class RootFactory(object):
    """Context factory which sets default acl"""
    __name__ = ''
    __acl__ = [
        (Allow, Authenticated, 'view'),
        (Deny, Everyone, ALL_PERMISSIONS)
    ]

    def __init__(self, request):
        self.request = request
        self.extjsurl()
        self.user()

    def extjsurl(self):
        """Adds extjsroot url as request.extjsroot
        Uses 'extjsroot' setting
        """
        extjsroot = 'magmaweb:static/'
        extjsroot = extjsroot + self.request.registry.settings['extjsroot']
        self.request.extjsroot = self.request.static_url(extjsroot)

    def user(self):
        """Adds :class:User as request.user

        Returns None when there is no authenticatied userid or
        if user is not found in user db
        """

        def userid2user(request):
            userid = unauthenticated_userid(request)
            if userid is not None:
                return DBSession().query(User).get(userid)
            else:
                return None

        # make request.user property lazy
        self.request.set_property(userid2user, 'user', reify=True)

class JobIdFactory(RootFactory):
    """Context factory which creates :class:magmaweb.job.Job
    as context of request
    """
    def __init__(self, request):
        super(JobIdFactory, self).__init__(request)
        self.job_factory = make_job_factory(request.registry.settings)

    def __getitem__(self, jobid):
        try:
            job = self.job_factory.fromId(uuid.UUID(jobid))
            # for acl inheritance add a parent,
            # this is not the parent job where this job is derived from
            job.__parent__ = self
            # owner may run calculations
            job.__acl__ = [(Allow, job.owner, 'run')]
        except JobNotFound:
            # TODO to fetch job state job's db can be absent
            raise HTTPNotFound()
        return job
