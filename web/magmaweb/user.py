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
from pyramid.security import Allow, Deny, Everyone, ALL_PERMISSIONS, Authenticated
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

    def __init__(self, userid, displayname, email):
        self.userid = userid
        self.displayname = displayname
        self.email = email


class JobMeta(Base):
    """Job metadata"""
    __tablename__ = 'job'
    jobid = Column(UUIDType, primary_key=True)
    description = Column(Unicode) # Kept in sync with Run.description
    parentjobid = Column(UUIDType, ForeignKey('job.jobid'))
    state = Column(Unicode) # queued/running/ready etc.
    owner = Column(Unicode, ForeignKey('user.userid'))
    children = relationship('JobMeta', backref=backref('parent', remote_side=[jobid]))

    def __init__(self, jobid, description='', parentjobid=None, state=None, owner=None):
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

    def extjsurl(self):
        """Adds extjsroot url to request using extjsroot setting"""
        extjsroot = 'magmaweb:static/'+self.request.registry.settings['extjsroot']
        self.request.extjsroot = self.request.static_url(extjsroot)

class JobIdFactory(RootFactory):
    """Context factory which creates DataSet as context of request"""
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

def get_jobs(owner):
    """Retrieve jobs belonging to owner"""
    jobs = []
    for jobmeta in DBSession().query(JobMeta).filter(JobMeta.owner==owner):
        jobs.append({'id':str(jobmeta.jobid), "description": jobmeta.description})
    return jobs
