import transaction

from sqlalchemy import Column
from sqlalchemy import Integer
from sqlalchemy import Unicode
from sqlalchemy import Float
from sqlalchemy import Boolean

from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker

from zope.sqlalchemy import ZopeTransactionExtension

DBSession = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))
Base = declarative_base()

class Metabolite(Base):
    __tablename__ = 'metabolites'
    id = Column(Integer, primary_key=True)
    mol = Column(Unicode, unique=True)
    level = Column(Integer)
    probability = Column(Float)
    reactionsequence = Column(Unicode)
    smiles = Column(Unicode)
    molformula = Column(Unicode)
    isquery = Column(Boolean)
    origin = Column(Integer)
    nhits = Column(Integer)

def initialize_sql(engine):
    DBSession.configure(bind=engine)
    Base.metadata.bind = engine
