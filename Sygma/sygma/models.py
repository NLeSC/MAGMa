import transaction

from sqlalchemy import Column
from sqlalchemy import Integer
from sqlalchemy import Unicode
from sqlalchemy import Float
from sqlalchemy import Boolean
from sqlalchemy import ForeignKey

from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import relationship, backref

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
    origin = Column(Unicode)
    nhits = Column(Integer)

class Scan(Base):
    __tablename__ = 'scans'
    scanid = Column(Integer, primary_key=True)
    mslevel = Column(Integer)
    rt = Column(Float)
    lowmz = Column(Float)
    highmz = Column(Float)
    basepeakmz = Column(Float)
    basepeakintensity = Column(Float)
    totioncurrent = Column(Integer)
    precursormz = Column(Float)
    precursorintensity = Column(Float)
    precursorscanid = Column(Integer, ForeignKey('scans.scanid'))
    peaks = relationship('Peak', backref='scan')
    precursor = relationship('Scan')

class Peak(Base):
    __tablename__ = 'peaks'
    scanid = Column(Integer, ForeignKey('scans.scanid'), primary_key=True)
    mz = Column(Float, primary_key=True)
    intensity = Column(Float)

def initialize_sql(engine):
    DBSession.configure(bind=engine)
    Base.metadata.bind = engine
