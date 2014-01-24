"""
Sqlalchemy models for magma result database
"""

import json
from sqlalchemy import Column
from sqlalchemy import Integer
from sqlalchemy import Unicode
from sqlalchemy import Float
from sqlalchemy import Boolean
from sqlalchemy import ForeignKey
from sqlalchemy import TypeDecorator

from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy.orm import relationship, backref
from sqlalchemy.schema import ForeignKeyConstraint

Base = declarative_base()

class ReactionSequence(TypeDecorator):
    """List of reactions.

    Reactions are grouped by relation to the current row or molecule:

    * products, reactions which have current row as reactant
    * reactants, reactions which have current row as product

    Reactions is a dict with key as the reaction name and the value a dict with keys:

    * nr, number of molecules which are product/reactant
        of current molecule with this reaction
    * nrp, number of molecules which are product/reactant
        of current molecule with this reaction and which have been matched to at least one scan

    Example:

    .. code-block:: javascript

          {
            'reactants': {
               'esterase': {'nr': 123, 'nrp': 45}
            },
            'products': {
               'sulfation': {'nr': 678, 'nrp': 90}
            }
          }

    Stored in database as json serialized string.
    """
    impl = Unicode

    def process_bind_param(self, value, dialect):
        if value is not None:
            value = json.dumps(value)
        return value

    def process_result_value(self, value, dialect):
        if value is u'':
            value = u'{}'
        if value is not None:
            try:
                value = json.loads(value)
            except ValueError:
                value = {}
        return value

class Metabolite(Base):
    """Metabolite model for metabolites table"""
    __tablename__ = 'metabolites'
    metid = Column(Integer, primary_key=True, autoincrement=True) #: Id of a metabolite
    mol = Column(Unicode) #: molfile as string
    level = Column(Integer)
    probability = Column(Float)
    # A newline seperated list of reactions
    reactionsequence = Column(ReactionSequence, default={})
    smiles = Column(Unicode, unique=True) #: Smiles string
    molformula = Column(Unicode) #: Molecular formula
    isquery = Column(Boolean) #: Whether metabolite was given as query or is a result a of reaction
    origin = Column(Unicode) #: Name of molecule
    nhits = Column(Integer)
    mim = Column(Float) #: Monoisotopic mass
    natoms = Column(Integer) #: Number of non-hydrogen atoms
    logp = Column(Float) #: Calculated logP
    reference = Column(Unicode)
    fragments = relationship('Fragment', backref='metabolite') #: each metabolite is fragmented into fragments

class Reaction(Base):
    """Reaction model for reactions table"""
    __tablename__ = 'reactions'
    reactid = Column(Integer, primary_key=True, autoincrement=True) #: Id of a reaction
    reactant = Column(Integer)
    product = Column(Integer)
    name = Column(Unicode)

def fill_molecules_reactions(session):
    """Fills the reactionsequence column in the molecules table with info from reactions table.

    The molecules query will become to complex when reactionsequence is queried from reactions table.
    So we fill reaction sequence with a json serialized struct which can be used during rendering.

    from magmaweb.job import JobFactory
    factory = JobFactory('data/jobs')
    session = factory._makeJobSession('58f05077-aad8-4fc9-a497-310495ab8b62')
    from magmaweb.models import fill_molecules_reactionsequence
    fill_molecules_reactionsequence(session)

    """
    from sqlalchemy.sql import func

    reactions = {}
    for metid, rname, nr in session.query(Reaction.product, Reaction.name, func.count('*')).group_by(Reaction.product, Reaction.name):
        if metid not in reactions:
            reactions[metid] = {}
        if 'reactants' not in reactions[metid]:
            reactions[metid]['reactants'] = {}
        reactions[metid]['reactants'][rname] = {'nr': nr}

    for metid, rname, nrp in session.query(Reaction.product, Reaction.name, func.count('*')).join(Metabolite, Metabolite.metid == Reaction.reactant).filter(Metabolite.nhits > 0).group_by(Reaction.product, Reaction.name):
        # dont need checks for keys because query above is always superset of this query
        reactions[metid]['reactants'][rname]['nrp'] = nrp

    for metid, rname, nr in session.query(Reaction.reactant, Reaction.name, func.count('*')).group_by(Reaction.reactant, Reaction.name):
        if metid not in reactions:
            reactions[metid] = {}
        if 'products' not in reactions[metid]:
            reactions[metid]['products'] = {}
        reactions[metid]['products'][rname] = {'nr': nr}

    for metid, rname, nrp in session.query(Reaction.reactant, Reaction.name, func.count('*')).join(Metabolite, Metabolite.metid == Reaction.product).filter(Metabolite.nhits > 0).group_by(Reaction.reactant, Reaction.name):
        # dont need checks for keys because query above is always superset of this query
        reactions[metid]['products'][rname]['nrp'] = nrp

    for mol in session.query(Metabolite):
        if mol.metid in reactions:
            reaction = reactions[mol.metid]
        else:
            reaction = {}
        mol.reactionsequence = reaction
    session.commit()

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
    precursorscanid = Column(Integer, ForeignKey('scans.scanid'),index=True) #: Parent scan identifier
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
    assigned_metid = Column(Integer, ForeignKey('metabolites.metid'),index=True) # which metabolite is assigned to this peak

class Fragment(Base):
    """Fragment model for fragments table"""
    __tablename__ = 'fragments'
    fragid = Column(Integer, primary_key=True, autoincrement=True) #: Fragment identifier
    metid = Column(Integer, ForeignKey('metabolites.metid'),index=True) #: Metabolite identifier
    scanid = Column(Integer, ForeignKey('scans.scanid'),index=True) #: Scan identifier
    mz = Column(Float) #: m/z of peak in scan
    mass = Column(Float) #: Mass of fragment in Dalton, corrected with h delta
    score = Column(Float) #: Score of how good the metabolite fragment matches the mass spectras
    parentfragid = Column(Integer, ForeignKey('fragments.fragid'),index=True)
    atoms = Column(Unicode) #: Atom indices of metabolite which are the fragment is a comma seperated list, starting with 0
    deltah = Column(Float)
    deltappm = Column(Float)
    inchikey = Column(Unicode)
    formula = Column(Unicode)
    #: A fragment can have child fragments
    children = relationship('Fragment', backref=backref('parent', remote_side=[fragid]), lazy='joined', join_depth=1)
    __table_args__ = (ForeignKeyConstraint(['scanid', 'mz'],
                                           ['peaks.scanid', 'peaks.mz']
                                           ),
                      {}
                      )

class Run(Base):
    """Run model for run table"""
    __tablename__ = 'run'
    runid = Column(Integer, primary_key=True, autoincrement=True) #: Run identifier
    description = Column(Unicode) #: Description of the run

    # SyGMa parameters, TODO remove: metabolism type info will be part of reacton sequence of metabolites
    n_reaction_steps = Column(Integer) #: Maximum number of reaction steps applied to reactants
    metabolism_types = Column(Unicode) #: Comma separated list of metabolism types, like "phase1"
    # ms data parsing parameters
    ms_filename = Column(Unicode)
    abs_peak_cutoff = Column(Float) #: abs intensity threshold for storing peaks in database
    max_ms_level = Column(Integer) #: maximum ms level to be included in the analysis

    # parameters for matching metabolites and fragments with peaks
    ionisation_mode = Column(Integer)
    skip_fragmentation = Column(Boolean)
    max_broken_bonds = Column(Integer) #: max number of bonds broken in substructures generated from metabolites
    max_water_losses = Column(Integer) #: max number of additional neutral water losses
    ms_intensity_cutoff = Column(Float) #: Absolute intensity minimum of lvl1 scan peaks which are matched with metabolites
    msms_intensity_cutoff = Column(Float) #: Ratio of basepeak intensity
    mz_precision = Column(Float) #: precision for matching a metabolite mim to m/z of a peak
    mz_precision_abs = Column(Float) #: precision for matching a metabolite mim to m/z of a peak
    precursor_mz_precision = Column(Float) #: precision for matching precursor mz with peak mz in parent scan
    use_all_peaks = Column(Boolean)
