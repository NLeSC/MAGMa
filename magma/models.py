"""
Sqlalchemy models for msygma result database
"""

from sqlalchemy import Column
from sqlalchemy import Integer
from sqlalchemy import Unicode
from sqlalchemy import Float
from sqlalchemy import Boolean
from sqlalchemy import ForeignKey

from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy.orm import relationship, backref

Base = declarative_base()

class Metabolite(Base):
    """Metabolite model for metabolites table"""
    __tablename__ = 'metabolites'
    metid = Column(Integer, primary_key=True, autoincrement=True) #: Id of a metabolite
    mol = Column(Unicode) #: molfile as string
    level = Column(Integer)
    probability = Column(Float)
    reactionsequence = Column(Unicode) #: A newline seperated list of reactions
    smiles = Column(Unicode, unique=True) #: Smile string
    molformula = Column(Unicode) #: Molecular formula
    isquery = Column(Boolean) #: Whether metabolite was given as query or is a result a of reaction
    origin = Column(Unicode) #: Name of molecule
    nhits = Column(Integer)
    mim = Column(Float) #: Monoisotopic mass
    fragments = relationship('Fragment', backref='metabolite') #: each metabolite is fragmented into fragments

class Scan(Base):
    """Scan model for scans table"""
    __tablename__ = 'scans'
    scanid = Column(Integer, primary_key=True) #: Id of a scan
    mslevel = Column(Integer) #: Level of ms, starts with 1
    rt = Column(Float) #: Retention time on which scan was taken
    lowmz = Column(Float) #: Lowest m/z
    highmz = Column(Float) #: Highest m/z
    basepeakmz = Column(Float) #: M/z with highest intensity
    basepeakintensity = Column(Float) #: Highest intensity
    totioncurrent = Column(Float)
    precursormz = Column(Float) #: m/z of precursor (which was fragmented resulting in this scan)
    precursorintensity = Column(Float) #: intensity belonging to precursorintensity
    precursorscanid = Column(Integer, ForeignKey('scans.scanid')) #: Parent scan identifier
    peaks = relationship('Peak', backref='scan') #: Each scan has many peaks
    #: A scan can have child product scans and a parent precursor scan
    products = relationship('Scan', backref=backref('precursor', remote_side=[scanid]))
    fragments = relationship('Fragment', backref='scan') #: Fragments can be found on a scan

class Peak(Base):
    """Peak model for peaks table"""
    __tablename__ = 'peaks'
    scanid = Column(Integer, ForeignKey('scans.scanid'), primary_key=True) #: Scan identifier to which peaks belongs
    mz = Column(Float, primary_key=True) #: m/z of peak (x-coordinate)
    intensity = Column(Float) #: Intensity of peak (y-coordinate)

class Fragment(Base):
    """Fragment model for fragments table"""
    __tablename__ = 'fragments'
    fragid = Column(Integer, primary_key=True, autoincrement=True) #: Fragment identifier
    metid = Column(Integer, ForeignKey('metabolites.metid')) #: Metabolite identifier
    scanid = Column(Integer, ForeignKey('scans.scanid'), ForeignKey('peaks.scanid')) #: Scan identifier
    mz = Column(Float, ForeignKey('peaks.scanid')) #: m/z of peak in scan
    mass = Column(Float) #: Mass of fragment in Dalton, corrected with h delta
    score = Column(Float) #: Score of how good the metabolite fragment matches the mass spectras
    parentfragid = Column(Integer, ForeignKey('fragments.fragid'))
    atoms = Column(Unicode) #: Atom indices of metabolite which are the fragment is a comma seperated list, starting with 0
    deltah = Column(Float)
    #: A fragment can have child fragments
    children = relationship('Fragment', backref=backref('parent', remote_side=[fragid]), lazy='joined', join_depth=1)

class Run(Base):
    """Run model for run table"""
    __tablename__ = 'run'
    # TODO run consists of 1 row so doesn't need a pk
    # but sqlalchemy requires one

    # SyGMa parameters, TODO remove: metabolism type info will be part of reacton sequence of metabolites
    n_reaction_steps = Column(Integer, primary_key=True) #: Maximum number of reaction steps applied to reactants
    metabolism_types = Column(Unicode) #: Comma separated list of metabolism types, like "phase1"

    # ms data parsing parameters
    ms_filename = Column(Unicode)
    abs_peak_cutoff = Column(Float) #: abs intensity threshold for storing peaks in database
    rel_peak_cutoff = Column(Float) #: fraction of basepeak intensity threshold for storing peaks in database
    max_ms_level = Column(Integer) #: maximum ms level to be included in the analysis

    # parameters for matching metabolites and fragments with peaks
    ionisation_mode = Column(Integer)
    use_fragmentation = Column(Boolean)
    max_broken_bonds = Column(Integer) #: max number of bonds broken in substructures generated from metabolites
    ms_intensity_cutoff = Column(Float) #: Absolute intensity minimum of lvl1 scan peaks which are matched with metabolites
    msms_intensity_cutoff = Column(Float) #: Ratio of basepeak intensity
    mz_precision = Column(Float) #: precision for matching a metabolite mim to m/z of a peak
    precursor_mz_precision = Column(Float) #: precision for matching precursor mz with peak mz in parent scan
    use_msms_only = Column(Boolean)
