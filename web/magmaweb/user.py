"""Module for authentication, authorization.

Contains db schema for users and jobs.
"""
import uuid
import datetime
import bcrypt
from sqlalchemy import Column
from sqlalchemy import Unicode
from sqlalchemy import String
from sqlalchemy import ForeignKey
from sqlalchemy import DateTime
from sqlalchemy import Boolean
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import backref
from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import validates
from sqlalchemy.types import TypeDecorator
from zope.sqlalchemy import ZopeTransactionExtension
from pyramid.httpexceptions import HTTPNotFound
from pyramid.security import Allow, Deny, Everyone, ALL_PERMISSIONS
from pyramid.security import Authenticated, authenticated_userid

from magmaweb.job import make_job_factory, JobNotFound

DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))

Base = declarative_base()


def init_user_db(engine, create=True, fill=True):
    """Binds engine with DBSession and Mappings.

    'engine' is a :class:`sqlalchemy.engine.base.Engine`.
    Set 'create' to False to skip createing tables.
    Set 'fill' to False to skip adding 'joblauncher' user.
    """
    DBSession.configure(bind=engine)
    Base.metadata.bind = engine
    if create:
        Base.metadata.create_all(engine)
    if fill:
        # TODO register all found jobs
        pass


class UUIDType(TypeDecorator):
    """ Stores UUID as Unicode string in database"""
    impl = Unicode

    def process_bind_param(self, value, dialect):
        if value is None:
            return None
        return unicode(value)

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
            self.password = password

    @validates('password')
    def _hash_password(self, key, clear_password):
        salt = bcrypt.gensalt(12)
        return bcrypt.hashpw(clear_password, salt)

    def validate_password(self, clear_password):
        """Validates unencrypted password 'clear_password' matches

        Returns True when matches or False when not
        """
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
        """Fetch :class:`User` by `userid`"""
        return DBSession().query(cls).get(userid)

    @classmethod
    def add(cls, user):
        """Adds :class:`User` to db"""
        session = DBSession()
        session.add(user)
        session.flush()


class JobMeta(Base):
    """Job metadata"""
    __tablename__ = 'job'
    jobid = Column(UUIDType, primary_key=True)
    description = Column(Unicode)  # Kept in sync with Run.description
    ms_filename = Column(Unicode)  # Kept in sync with Run.ms_filename
    parentjobid = Column(UUIDType, ForeignKey('job.jobid'))
    state = Column(Unicode, nullable=False)  # queued/running/ready etc.
    owner = Column(Unicode, ForeignKey('user.userid'), nullable=False)
    created_at = Column(DateTime, nullable=False)
    is_public = Column(Boolean, default=False)
    children = relationship('JobMeta',
                            backref=backref('parent', remote_side=[jobid]))

    def __init__(self, jobid, owner,
                 description=u'', parentjobid=None,
                 state=u'STOPPED', ms_filename='',
                 created_at=None,
                 is_public=False):
        self.jobid = jobid
        self.owner = owner
        self.description = description
        self.parentjobid = parentjobid
        self.state = state
        self.ms_filename = ms_filename
        self.is_public = is_public
        if created_at is None:
            self.created_at = datetime.datetime.utcnow()
        else:
            self.created_at = created_at

    @classmethod
    def by_id(cls, jobid):
        """Fetch :class:`JobMeta` by `jobid`"""
        return DBSession().query(cls).get(jobid)

    @classmethod
    def add(cls, jobmeta):
        """Adds :class:JobMeta to db"""
        session = DBSession()
        session.add(jobmeta)
        # commit is normally done at the end of request handling
        # Exceptions cause transactions to be aborted
        # so force commit here
        session.flush()

    @classmethod
    def delete(cls, jobmeta):
        """Deletes :class:JobMeta from db"""
        session = DBSession()
        session.delete(jobmeta)


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
        """Adds :class:`User` as request.user

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
    """Context factory which creates :class:`magmaweb.job.Job`
    as context of request
    """
    def __init__(self, request):
        super(JobIdFactory, self).__init__(request)
        self.job_factory = make_job_factory(request.registry.settings)

    def __getitem__(self, jobid):
        try:
            job = self.job_factory.fromId(uuid.UUID(jobid))
            # owner may view and run calculations or perform changes
            monitor_user = self.request.registry.settings['monitor_user']
            job.__acl__ = [(Allow, job.owner, ('run', 'view')),
                           # monitor user may monitor this job
                           (Allow, monitor_user, 'monitor'),
                           (Deny, Everyone, ALL_PERMISSIONS),
                           ]
            # if job has been marked public it can be viewed
            # by all authenticated users if they know the job url.
            if job.is_public:
                job.__acl__.insert(0, (Allow, Authenticated, 'view'))
        except JobNotFound:
            # TODO to fetch job state job's db can be absent
            raise HTTPNotFound()
        return job
