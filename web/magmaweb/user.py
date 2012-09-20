from sqlalchemy import Column
from sqlalchemy import Unicode
from sqlalchemy import ForeignKey
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
    roles = relationship('UserRole',backref='user')

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
    children = relationship('Job', backref=backref('parent', remote_side=[jobid]))

    def __init__(self, jobid, description='', parentjobid=None, state=None):
        self.jobid = jobid
        self.description = description
        self.parentjobid = parentjobid
        self.state = state

class JobUser(Base):
    """Which user has which roles on a job"""
    __tablename__ = 'job_user'
    userid = Column(Unicode, ForeignKey('user.userid'), primary_key=True)
    jobid = Column(Unicode, ForeignKey('job.jobid'), primary_key=True)
    role = Column(Unicode)

    def __init__(self, userid, jobid, role):
        self.userid = userid
        self.jobid = jobid
        self.role = role

class UserRole(Base):
    """Roles of a user"""
    __tablename__ = 'user_role'
    userid = Column(Unicode, ForeignKey('user.userid'), primary_key=True)
    role = Column(Unicode, primary_key=True)

    def __init__(self, userid, role):
        self.userid = userid
        self.role = role

def groupfinder(userid, request):
    roles = DBSession.query(UserRole).filter(UserRole.userid==userid).all()
    return [r.role for r in roles]

class RootFactory(object):
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
    def __init__(self, request):
        super(JobIdFactory, self).__init__(request)
        self.job_factory = make_job_factory(request.registry.settings)

    def __getitem__(self, jobid):
        try:
            job = self.job_factory.fromId(jobid)
        except JobNotFound:
            raise HTTPNotFound()
        owner = DBSession.query(JobUser).filter(JobUser.jobid==jobid).filter(JobUser.role=='owner').one()
        job.__acl__ = [
            (Allow, owner.userid, 'run')
        ]
        job.__parent__ = self
        #TODO: add users/groups with other permission on job to acl
        return job

def set_job_owner(jobid, userid):
    """ userid=unauthenticated_userid(self.request)"""
    # TODO: Only one user should be owner, revoke others ownership
    grant('owner', jobid, userid)

def grant(role, jobid, userid):
    """ GRANT role ON job TO who"""
    ju = JobUser(userid, jobid, role)
    DBSession().add(ju)

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
