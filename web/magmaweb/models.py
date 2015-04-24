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
from sqlalchemy.sql import func

Base = declarative_base()


class ReactionSequence(TypeDecorator):
    """List of reactions.

    Reactions are grouped by relation to the current row or molecule:

    * products, reactions which have current row as reactant
    * reactants, reactions which have current row as product

    Reactions is a dict with key as the reaction name and
    the value a dict with keys:

    * nr, number of molecules which are product/reactant
        of current molecule with this reaction
    * nrp, number of molecules which are product/reactant
        of current molecule with this reaction and
        which have been matched to at least one scan

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
            value = unicode(json.dumps(value))
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


class Molecule(Base):
    """Molecule model for molecules table"""
    __tablename__ = 'molecules'
    # Id of a molecule
    molid = Column(Integer, primary_key=True, autoincrement=True)
    # molfile as string
    mol = Column(Unicode)
    # A newline seperated list of reactions
    reactionsequence = Column(ReactionSequence, default={})
    # Inchikey first 14 chars
    inchikey14 = Column(Unicode, unique=True)
    # Smile string
    smiles = Column(Unicode)
    # Molecular formula
    formula = Column(Unicode)
    # Whether molecule was given as query or is a result a of reaction
    predicted = Column(Boolean)
    # Name of molecule
    name = Column(Unicode)
    # Number of lvl1 scans fragments are found for this molecule
    nhits = Column(Integer)
    # Monoisotopic mass
    mim = Column(Float)
    # Number of non-hydrogen atoms
    natoms = Column(Integer)
    # Calculated logP
    logp = Column(Float)
    refscore = Column(Float)
    reference = Column(Unicode)
    # each molecule is fragmented into fragments
    fragments = relationship('Fragment', backref='molecule')


class Reaction(Base):
    """Reaction model for reactions table"""
    __tablename__ = 'reactions'
    # Id of a reaction
    reactid = Column(Integer, primary_key=True, autoincrement=True)
    reactant = Column(Integer, ForeignKey('molecules.molid'))
    product = Column(Integer, ForeignKey('molecules.molid'))
    name = Column(Unicode)


def fill_molecules_reactions_reactants(session, reactions):
    qrnr = session.query(Reaction.product, Reaction.name, func.count('*'))
    qrnr = qrnr.group_by(Reaction.product, Reaction.name)
    for molid, rname, nr in qrnr:
        if molid not in reactions:
            reactions[molid] = {}
        if 'reactants' not in reactions[molid]:
            reactions[molid]['reactants'] = {}
        reactions[molid]['reactants'][rname] = {'nr': nr}

    qrnrp = session.query(Reaction.product, Reaction.name, func.count('*'))
    qrnrp = qrnrp.join(Molecule, Molecule.molid == Reaction.reactant)
    qrnrp = qrnrp.filter(Molecule.nhits > 0)
    qrnrp = qrnrp.group_by(Reaction.product, Reaction.name)
    for molid, rname, nrp in qrnrp:
        # dont need checks for keys because query above is always superset of
        # this query
        reactions[molid]['reactants'][rname]['nrp'] = nrp


def fill_molecules_reactions_products(session, reactions):
    qpnr = session.query(Reaction.reactant, Reaction.name, func.count('*'))
    qpnr = qpnr.group_by(Reaction.reactant, Reaction.name)
    for molid, rname, nr in qpnr:
        if molid not in reactions:
            reactions[molid] = {}
        if 'products' not in reactions[molid]:
            reactions[molid]['products'] = {}
        reactions[molid]['products'][rname] = {'nr': nr}

    qpnrp = session.query(Reaction.reactant, Reaction.name, func.count('*'))
    qpnrp = qpnrp.join(Molecule, Molecule.molid == Reaction.product)
    qpnrp = qpnrp.filter(Molecule.nhits > 0)
    qpnrp = qpnrp.group_by(Reaction.reactant, Reaction.name)
    for molid, rname, nrp in qpnrp:
        # dont need checks for keys because query above is always superset of
        # this query
        reactions[molid]['products'][rname]['nrp'] = nrp


def fill_molecules_reactions(session):
    """Fills the reactionsequence column in the molecules table with
    info from reactions table.

    The molecules query will become to complex when
    reactionsequence is queried from reactions table.
    So we fill reaction sequence with a json serialized struct
    which can be used during rendering.

    from magmaweb.job import JobFactory
    factory = JobFactory('data/jobs')
    session = factory._makeJobSession('58f05077-aad8-4fc9-a497-310495ab8b62')
    from magmaweb.models import fill_molecules_reactionsequence
    fill_molecules_reactionsequence(session)

    """

    reactions = {}
    fill_molecules_reactions_reactants(session, reactions)
    fill_molecules_reactions_products(session, reactions)

    for mol in session.query(Molecule):
        if mol.molid in reactions:
            reaction = reactions[mol.molid]
        else:
            reaction = {}
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
    # A scan can have child product scans and a parent precursor scan
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
    # which molecule is assigned to this peak
    assigned_molid = Column(Integer, ForeignKey('molecules.molid'), index=True)


class Fragment(Base):
    """Fragment model for fragments table"""
    __tablename__ = 'fragments'
    # Fragment identifier
    fragid = Column(Integer, primary_key=True, autoincrement=True)
    # Molecule identifier
    molid = Column(Integer, ForeignKey('molecules.molid'), index=True)
    # Scan identifier
    scanid = Column(Integer, ForeignKey('scans.scanid'), index=True)
    # m/z of peak in scan
    mz = Column(Float)
    # Mass of fragment in Dalton, corrected with h delta
    mass = Column(Float)
    # Score of how good the molecule fragment matches the mass spectras
    score = Column(Float)
    # From which fragment this fragment is a fragment of
    parentfragid = Column(Integer, ForeignKey('fragments.fragid'))
    # Atom indices of molecule which are the fragment
    # is a comma seperated list, starting with 0
    atoms = Column(Unicode)
    deltah = Column(Float)
    # (mz+deltah*1.007825032-mass)/(mz*1e6)  as deltappm
    deltappm = Column(Float)
    smiles = Column(Unicode)
    # molecular formula of fragment
    formula = Column(Unicode)
    # A fragment can have child fragments
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

    # ms data parsing parameters
    ms_filename = Column(Unicode)
    # abs intensity threshold for storing peaks in database
    abs_peak_cutoff = Column(Float)
    # maximum ms level to be included in the analysis
    max_ms_level = Column(Integer)

    # parameters for matching molecules and fragments with peaks
    ionisation_mode = Column(Integer)
    skip_fragmentation = Column(Boolean)
    # max number of bonds broken in substructures generated from molecules
    max_broken_bonds = Column(Integer)
    # max number of additional neutral water losses
    max_water_losses = Column(Integer)
    # Absolute intensity minimum of lvl1 scan peaks
    # which are matched with molecules
    ms_intensity_cutoff = Column(Float)
    # Ratio of basepeak intensity
    msms_intensity_cutoff = Column(Float)
    # precision for matching a molecule mim to m/z of a peak (in ppm)
    mz_precision = Column(Float)
    # precision for matching a molecule mim to m/z of a peak (in Da)
    mz_precision_abs = Column(Float)
    # precision for matching precursor mz with peak mz in parent scan
    precursor_mz_precision = Column(Float)
    use_all_peaks = Column(Boolean)
