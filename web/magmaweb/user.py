import uuid
import bcrypt
from sqlalchemy import Column
from sqlalchemy import Unicode
from sqlalchemy import String
from sqlalchemy import ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker
from sqlalchemy.types import TypeDecorator
from zope.sqlalchemy import ZopeTransactionExtension #@UnresolvedImport
from pyramid.httpexceptions import HTTPNotFound
from pyramid.security import Allow, Deny, Everyone, ALL_PERMISSIONS
from pyramid.security import Authenticated, authenticated_userid

from magmaweb.job import make_job_factory, JobNotFound

DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))

Base = declarative_base()


class UUIDType(TypeDecorator):
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
    __tablename__ = 'user'
    userid = Column(Unicode, primary_key=True)
    displayname = Column(Unicode)
    email = Column(Unicode)
    password = Column(String)  # bcrypted hashed pw
    jobs = relationship("JobMeta", order_by="JobMeta.jobid")

    def __init__(self, userid, displayname, email, password=None):
        self.userid = userid
        self.displayname = displayname
        self.email = email
        if password is not None:
            self._set_password(password)

    def _set_password(self, clear_password):
        salt = bcrypt.gensalt(12)
        self.password = bcrypt.hashpw(clear_password, salt)

    def validate_password(self, clear_password):
        return bcrypt.hashpw(clear_password, self.password) == self.password

    def __repr__(self):
        return '<User({!r}, {!r}, {!r})>'.format(self.userid,
                                                self.displayname,
                                                self.email)

    def __eq__(self, other):
        return (self.userid == other.userid and
            self.displayname == other.displayname and
            self.email == other.email and
            self.password == other.password)

    @classmethod
    def by_id(cls, userid):
        return DBSession().query(cls).get(userid)


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
            userid = authenticated_userid(request)
            if userid is not None:
                return User.by_id(userid)
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
