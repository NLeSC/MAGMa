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

    Reactions are grouped by relation to of the current row or molecule:

    * reactantof, reactions which have current row as reactant
    * productof, reactions which have current as product

    Reactions is a dict with key as the reaction name and the value a dict with keys:

    * nr, number of molecules which are product/reactant
        of current molecule with this reaction
    * nrp, number of molecules which are product/reactant
        of current molecule with this reaction and which have been matched to at least one scan

    Example:

          {
            'reactantof': [{
               'esterase': {'nr': 123, 'nrp': 45}
            }],
            'productof': [{
               'theogallin': {'nr': 678, 'nrp': 90}
            }]
          }

    Stored in database as json serialized string.
    """
    impl = Unicode

    def process_bind_param(self, value, dialect):
        if value is not None:
            value = json.dumps(value)

        return value

    def process_result_value(self, value, dialect):
        if value is '':
            value = '{}'
        if value is not None:
            value = json.loads(value)
        return value


class Metabolite(Base):
    """Metabolite model for metabolites table"""
    __tablename__ = 'metabolites'
    # Id of a metabolite
    metid = Column(Integer, primary_key=True, autoincrement=True)
    # molfile as string
    mol = Column(Unicode)
    level = Column(Integer)
    probability = Column(Float)
    # A newline seperated list of reactions
    reactionsequence = Column(ReactionSequence, default={})
    # Smile string
    smiles = Column(Unicode, unique=True)
    # Molecular formula
    molformula = Column(Unicode)
    # Whether metabolite was given as query or is a result a of reaction
    isquery = Column(Boolean)
    # Name of molecule
    origin = Column(Unicode)
    # Number of lvl1 scans fragments are found for this metabolite
    nhits = Column(Integer)
    # Monoisotopic mass
    mim = Column(Float)
    # Number of non-hydrogen atoms
    natoms = Column(Integer)
    # Calculated logP
    logp = Column(Float)
    reference = Column(Unicode)
    # each metabolite is fragmented into fragments
    fragments = relationship('Fragment', backref='metabolite')


class Reaction(Base):
    """Reaction model for reactions table"""
    __tablename__ = 'reactions'
    # Id of a reaction
    reactid = Column(Integer, primary_key=True, autoincrement=True)
    reactant = Column(Integer, ForeignKey('metabolites.metid'))
    product = Column(Integer, ForeignKey('metabolites.metid'))
    name = Column(Unicode)


def fill_molecules_reactionsequence(session):
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
        if 'productof' not in reactions[metid]:
            reactions[metid]['productof'] = {}
        if rname not in reactions[metid]['productof']:
            reactions[metid]['productof'][rname] = {}
        reactions[metid]['productof'][rname]['nr'] = nr

    for metid, rname, nrp in session.query(Reaction.product, Reaction.name, func.count('*')).join(Metabolite, Metabolite.metid == Reaction.reactant).filter(Metabolite.nhits > 0).group_by(Reaction.product, Reaction.name):
        # dont need checks for keys because query above is always superset of this query
        reactions[metid]['productof'][rname]['nrp'] = nrp

    for metid, rname, nr in session.query(Reaction.reactant, Reaction.name, func.count('*')).group_by(Reaction.reactant, Reaction.name):
        if metid not in reactions:
            reactions[metid] = {}
        if 'reactantof' not in reactions[metid]:
            reactions[metid]['reactantof'] = {}
        if rname not in reactions[metid]['reactantof']:
            reactions[metid]['reactantof'][rname] = {}
        reactions[metid]['reactantof'][rname]['nr'] = nr

    for metid, rname, nrp in session.query(Reaction.reactant, Reaction.name, func.count('*')).join(Metabolite, Metabolite.metid == Reaction.product).filter(Metabolite.nhits > 0).group_by(Reaction.reactant, Reaction.name):
        # dont need checks for keys because query above is always superset of this query
        reactions[metid]['reactantof'][rname]['nrp'] = nrp

    for mol in session.query(Metabolite):
        if mol.metid in reactions:
            reaction = reactions[mol.metid]
        else:
            reaction = {}
        if 'reactantof' in reaction:
            for reaction_name, counts in reaction['reactantof'].iteritems():
                if 'nrp' not in counts:
                    reaction['reactantof'][reaction_name]['nrp'] = 0
        if 'productof' in reaction:
            for reaction_name, counts in reaction['productof'].iteritems():
                if 'nrp' not in counts:
                    reaction['productof'][reaction_name]['nrp'] = 0
        mol.reactionsequence = reaction

    session.commit()


class Scan(Base):
    """Scan model for scans table"""
    __tablename__ = 'scans'
    # Id of a scan
    scanid = Column(Integer, primary_key=True)
    # Level of ms, starts with 1
    mslevel = Column(Integer)
    # Retention time on which scan was taken
    rt = Column(Float)
    #: Lowest m/z
    lowmz = Column(Float)
    #: Highest m/z
    highmz = Column(Float)
    #: M/z with highest intensity
    basepeakmz = Column(Float)
    #: Highest intensity
    basepeakintensity = Column(Float)
    totioncurrent = Column(Float)
    # m/z of precursor (which was fragmented resulting in this scan)
    precursormz = Column(Float)
    #: intensity belonging to precursormz
    precursorintensity = Column(Float)
    # Parent scan identifier
    precursorscanid = Column(Integer, ForeignKey('scans.scanid'))
    # Each scan has many peaks
    peaks = relationship('Peak', backref='scan')
    #: A scan can have child product scans and a parent precursor scan
    products = relationship('Scan', backref=backref('precursor',
                                                    remote_side=[scanid]))
    # Fragments can be found on a scan
    fragments = relationship('Fragment', backref='scan')


class Peak(Base):
    """Peak model for peaks table"""
    __tablename__ = 'peaks'
    # Scan identifier to which peaks belongs
    scanid = Column(Integer, ForeignKey('scans.scanid'), primary_key=True)
    # m/z of peak (x-coordinate)
    mz = Column(Float, primary_key=True)
    # Intensity of peak (y-coordinate)
    intensity = Column(Float)
    # which metabolite is assigned to this peak
    assigned_metid = Column(Integer, ForeignKey('metabolites.metid'))


class Fragment(Base):
    """Fragment model for fragments table"""
    __tablename__ = 'fragments'
    # Fragment identifier
    fragid = Column(Integer, primary_key=True, autoincrement=True)
    # Metabolite identifier
    metid = Column(Integer, ForeignKey('metabolites.metid'))
    # Scan identifier
    scanid = Column(Integer, ForeignKey('scans.scanid'))
    # m/z of peak in scan
    mz = Column(Float)
    # Mass of fragment in Dalton, corrected with h delta
    mass = Column(Float)
    # Score of how good the metabolite fragment matches the mass spectras
    score = Column(Float)
    # From which fragment this fragment is a fragment of
    parentfragid = Column(Integer, ForeignKey('fragments.fragid'))
    # Atom indices of metabolite which are the fragment
    # is a comma seperated list, starting with 0
    atoms = Column(Unicode)
    deltah = Column(Float)
    # (mz+deltah*1.007825032-mass)/(mz*1e6)  as deltappm
    deltappm = Column(Float)
    inchikey = Column(Unicode)
    # molecular formula of fragment
    formula = Column(Unicode)
    #: A fragment can have child fragments
    children_backref = backref('parent', remote_side=[fragid])
    children = relationship('Fragment', backref=children_backref,
                            lazy='joined', join_depth=1)
    __table_args__ = (ForeignKeyConstraint(['scanid', 'mz'],
                                           ['peaks.scanid', 'peaks.mz']
                                           ),
                      {}
                      )


class Run(Base):
    """Run model for run table"""
    __tablename__ = 'run'
    # Run identifier
    runid = Column(Integer, primary_key=True, autoincrement=True)
    # Description of the run
    description = Column(Unicode)

    # metabolize parameters

    # Maximum number of reaction steps applied to reactants
    n_reaction_steps = Column(Integer)

    # Comma separated list of metabolism types, like "phase1"
    metabolism_types = Column(Unicode)
    # TODO: remove metabolism type column
    # will be part of reacton sequence of metabolites

    # ms data parsing parameters
    ms_filename = Column(Unicode)
    # abs intensity threshold for storing peaks in database
    abs_peak_cutoff = Column(Float)
    # maximum ms level to be included in the analysis
    max_ms_level = Column(Integer)

    # parameters for matching metabolites and fragments with peaks
    ionisation_mode = Column(Integer)
    skip_fragmentation = Column(Boolean)
    # max number of bonds broken in substructures generated from metabolites
    max_broken_bonds = Column(Integer)
    # max number of additional neutral water losses
    max_water_losses = Column(Integer)
    # Absolute intensity minimum of lvl1 scan peaks
    # which are matched with metabolites
    ms_intensity_cutoff = Column(Float)
    # Ratio of basepeak intensity
    msms_intensity_cutoff = Column(Float)
    # precision for matching a metabolite mim to m/z of a peak (in ppm)
    mz_precision = Column(Float)
    # precision for matching a metabolite mim to m/z of a peak (in Da)
    mz_precision_abs = Column(Float)
    # precision for matching precursor mz with peak mz in parent scan
    precursor_mz_precision = Column(Float)
    use_all_peaks = Column(Boolean)
