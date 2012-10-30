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
from pyramid.httpexceptions import HTTPNotFound
from pyramid.security import Allow, Deny, Everyone, ALL_PERMISSIONS, Authenticated
from zope.sqlalchemy import ZopeTransactionExtension #@UnresolvedImport
from magmaweb.job import make_job_factory, JobNotFound

DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))

Base = declarative_base()

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


class Job(Base):
    """Job list for acl"""
    __tablename__ = 'job'
    jobid = Column(Unicode, primary_key=True) # Use uuid.int
    description = Column(Unicode) # Kept in sync with Run.description
    parentjobid = Column(Unicode, ForeignKey('job.jobid'))
    state = Column(Unicode) # queued/running/ready etc.
    owner = Column(Unicode, ForeignKey('user.userid'))
    children = relationship('Job', backref=backref('parent', remote_side=[jobid]))

    def __init__(self, jobid, description='', parentjobid=None, state=None, owner=None):
        self.jobid = jobid
        self.description = description
        self.parentjobid = parentjobid
        self.state = state
        self.owner = owner


class RootFactory(object):
    """Context factory which sets default acl"""
    __acl__ = [
        (Allow, Authenticated, 'view'),
        (Deny, Everyone, ALL_PERMISSIONS)
    ]

    def __init__(self, request):
        self.request = request
        self.extjsurl()

    def extjsurl(self):
        """Adds extjsroot url to request using extjsroot setting"""
        self.request.extjsroot = self.request.static_url('magmaweb:static/'+self.request.registry.settings['extjsroot'])

class JobIdFactory(RootFactory):
    """Context factory which creates Job as context of request"""
    def __init__(self, request):
        super(JobIdFactory, self).__init__(request)
        self.job_factory = make_job_factory(request.registry.settings)

    def __getitem__(self, jobid):
        try:
            job = self.job_factory.fromId(jobid)
        except JobNotFound:
            raise HTTPNotFound()
        ujob = DBSession.query(Job).filter(Job.jobid==jobid).one()
        job.__acl__ = [(Allow, ujob.owner, 'run')]
        job.__parent__ = self
        return job

def set_job_owner(jobid, userid):
    """ userid=unauthenticated_userid(self.request)"""
    session = DBSession()
    job = session.query(Job).get(jobid)
    job.owner = userid
    session.add(job)

def set_job_parent(jobid, parentjobid):
    """ Set Job.parentjobid with jobid==jobid"""
    session = DBSession()
    job = session.query(Job).get(jobid)
    job.parentjobid = parentjobid
    session.add(job)

def set_job_description(jobid, description):
    """ Set Job.description in user db"""
    session = DBSession()
    job = session.query(Job).get(jobid)
    job.description = description
    session.add(job)

def add_job(jobid):
    job = Job(jobid)
    DBSession().add(job)

def get_jobs(owner):
    """Retrieve jobs belonging to owner"""
    jobs = []
    for id,desc in DBSession().query(Job.jobid, Job.description).filter(Job.owner==owner).all():
        jobs.append({'id':id, "description": desc})
    return jobs
