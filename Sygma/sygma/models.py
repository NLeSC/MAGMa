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
    metid = Column(Integer, primary_key=True)
    mol = Column(Unicode, unique=True)
    level = Column(Integer)
    probability = Column(Float)
    reactionsequence = Column(Unicode)
    smiles = Column(Unicode)
    molformula = Column(Unicode)
    isquery = Column(Boolean)
    origin = Column(Unicode)
    nhits = Column(Integer)
    fragments = relationship('Fragment', backref='metabolite')

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
    fragments = relationship('Fragment', backref='scan')

class Peak(Base):
    __tablename__ = 'peaks'
    scanid = Column(Integer, ForeignKey('scans.scanid'), primary_key=True)
    mz = Column(Float, primary_key=True)
    intensity = Column(Float)

class Fragment(Base):
    __tablename__ = 'fragments'
    fragid = Column(Integer, primary_key=True)
    metid = Column(Integer, ForeignKey('metabolites.metid'))
    scanid = Column(Integer, ForeignKey('scans.scanid'), ForeignKey('peaks.scanid'))
    mz = Column(Float, ForeignKey('peaks.scanid'))
    mass = Column(Float)
    score = Column(Float)
    parentfragid = Column(Integer, ForeignKey('fragments.fragid'))
    atoms = Column(Unicode) # , seperated, starting with 0
    deltah = Column(Float)
    children = relationship('Fragment', backref=backref('parent', remote_side=[fragid]), lazy='joined', join_depth=1)

class Run(Base):
    __tablename__ = 'run'
    # TODO run consists of 1 row so doesn't need a pk
    # but sqlalchemy requires one
    n_reaction_steps = Column(Integer, primary_key=True)
    use_phase1 = Column(Boolean)
    use_phase2 = Column(Boolean)
    ms_filename = Column(Unicode)
    ionisation = Column(Unicode)
    use_fragmentation = Column(Boolean)
    ms_intensity_cutoff = Column(Float)
    msms_intensity_cutoff = Column(Float)
    mz_precision = Column(Float)
    use_msms_only = Column(Boolean)

def initialize_sql(engine):
    DBSession.configure(bind=engine)
    Base.metadata.bind = engine
