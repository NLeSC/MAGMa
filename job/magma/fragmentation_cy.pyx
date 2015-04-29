#### cython: profile=True
import numpy
cimport numpy
import pars
import os
from rdkit import Chem

ctypedef struct bonded_atom:
    int nbonds
    int[8] atoms

ctypedef struct bond_breaks_score_pair:
    int breaks
    float score

cdef class FragmentEngine(object):

    cdef unsigned long long new_fragment, template_fragment
    cdef int max_broken_bonds, max_water_losses, ionisation_mode, molcharge
    cdef bonded_atom[64] bonded_atoms
    cdef double[64] atom_masses
    cdef list neutral_loss_atoms
    cdef int nbonds, natoms, accept, skip_fragmentation
    cdef unsigned long long[128] bonds
    cdef float[128] bondscore
    cdef numpy.ndarray fragment_masses_np
    cdef list fragment_masses, fragment_info
    cdef int[64] atomHs
    cdef dict atom_elements
    cdef char * mol

    def __init__(self, mol, max_broken_bonds, max_water_losses, ionisation_mode, skip_fragmentation, molcharge):
        cdef unsigned long long bondbits
        cdef float bondscore
        cdef int x, a1, a2

        self.mol = mol
        try:
            mol = Chem.MolFromMolBlock(str(self.mol))
            self.accept = 1
            self.natoms = mol.GetNumAtoms()
        except:
            self.accept = 0
            return
        if self.natoms > 64:
            self.accept = 0
            return
        self.max_broken_bonds = max_broken_bonds
        self.max_water_losses = max_water_losses
        self.ionisation_mode = ionisation_mode
        self.skip_fragmentation = skip_fragmentation
        self.molcharge = molcharge
        self.nbonds = mol.GetNumBonds()
        self.neutral_loss_atoms = []
        self.atom_elements = {}
        self.new_fragment = 0
        self.template_fragment = 0
        self.fragment_masses = (
            (max_broken_bonds + max_water_losses) * 2 + 1) * [0]
        self.fragment_info = [[0, 0, 0]]

        for x in range(self.natoms):
            self.bonded_atoms[x].nbonds = 0
            atom = mol.GetAtomWithIdx(x)
            self.atomHs[x] = atom.GetNumImplicitHs() + atom.GetNumExplicitHs()
            self.atom_masses[x] = pars.mims[atom.GetSymbol()] + pars.Hmass * self.atomHs[x]
            self.atom_elements[x] = atom.GetSymbol()
            if self.atom_elements[x] == 'O' and self.atomHs[x] == 1 and len(atom.GetBonds()) == 1:
                self.neutral_loss_atoms.append(x)
            if self.atom_elements[x] == 'N' and self.atomHs[x] == 2 and len(atom.GetBonds()) == 1:
                self.neutral_loss_atoms.append(x)
        for x in range(self.nbonds):
            bond = mol.GetBondWithIdx(x)
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            self.bonded_atoms[a1].atoms[self.bonded_atoms[a1].nbonds] = a2
            self.bonded_atoms[a1].nbonds += 1
            self.bonded_atoms[a2].atoms[self.bonded_atoms[a2].nbonds] = a1
            self.bonded_atoms[a2].nbonds += 1
            bondbits = (1ULL << a1) | (1ULL << a2)
            bondscore = pars.typew[bond.GetBondType()] * \
                        pars.heterow[bond.GetBeginAtom().GetSymbol() != 'C' or bond.GetEndAtom().GetSymbol() != 'C']
            self.bonds[x] = bondbits
            self.bondscore[x] = bondscore

    cdef void extend(self, int atom):
        cdef int a, bonded_a
        cdef unsigned long long atombit
        for a in range(self.bonded_atoms[atom].nbonds):
            bonded_a = self.bonded_atoms[atom].atoms[a]
            atombit = 1ULL << bonded_a
            if atombit & self.template_fragment and not atombit & self.new_fragment:
                self.new_fragment = self.new_fragment | atombit
                self.extend(bonded_a)

    def generate_fragments(self):
        cdef unsigned long long fragment, frag
        cdef int atom, a, bonded_a
        cdef bond_breaks_score_pair bbsp
        cdef set all_fragments, total_fragments, current_fragments, new_fragments
        frag = (1ULL << self.natoms) - 1
        if self.natoms == 64:
            frag = 18446744073709551615  # =(1<<64)-1
        all_fragments = set([frag])
        total_fragments = set([frag])
        current_fragments = set([frag])
        new_fragments = set([frag])
        self.add_fragment(frag, self.calc_fragment_mass(frag), 0, 0)

        if self.skip_fragmentation:
            self.convert_fragments_table()
            return len(self.fragment_info)

        # generate fragments for max_broken_bond steps
        for step in range(self.max_broken_bonds):
            # loop of all fragments to be fragmented
            for fragment in current_fragments:
                # loop over all atoms
                for atom in range(self.natoms):
                    # in the fragment    
                    if (1ULL << atom) & fragment:
                        # remove the atom
                        self.template_fragment = fragment ^ (1ULL << atom)
                        list_ext_atoms = set([])
                        extended_fragments = set([])
                        # find all its neighbor atoms
                        for a in range(self.bonded_atoms[atom].nbonds):
                            bonded_a = self.bonded_atoms[atom].atoms[a]
                            # present in the fragment
                            if (1ULL << bonded_a) & self.template_fragment:
                                list_ext_atoms.add(bonded_a)
                        # in case of one bonded atom, the new fragment is the remainder of the old fragment
                        if len(list_ext_atoms) == 1:
                            extended_fragments.add(self.template_fragment)
                        else:
                           # otherwise extend each neighbor atom to a complete fragment
                            for a in list_ext_atoms:
                                # except when deleted atom is in a ring and a previous extended
                                # fragment already contains this neighbor atom, then
                                # calculate fragment only once
                                for frag in extended_fragments:
                                    if (1ULL << a) & frag:
                                        break
                                else:
                                    # extend atom to complete fragment
                                    self.new_fragment = 1ULL << a
                                    self.extend(a)
                                    extended_fragments.add(self.new_fragment)
                        for frag in extended_fragments:
                            # add extended fragments, if not yet present, to the collection
                            if frag not in all_fragments:
                                all_fragments.add(frag)
                                bbsp = self.score_fragment(frag)
                                if bbsp.breaks <= self.max_broken_bonds and bbsp.score < (pars.missingfragmentpenalty + 5):
                                    new_fragments.add(frag)
                                    total_fragments.add(frag)
                                    self.add_fragment(
                                        frag, self.calc_fragment_mass(frag), bbsp.score, bbsp.breaks)
            current_fragments = new_fragments
            new_fragments = set([])
        # number of OH losses
        for step in range(self.max_water_losses):
            # loop of all fragments
            for fi in self.fragment_info:
                # on which to apply neutral loss rules
                if fi[2] == self.max_broken_bonds + step:
                    fragment = fi[0]
                    # loop over all atoms in the fragment
                    for atom in self.neutral_loss_atoms:
                        if (1ULL << atom) & fragment:
                            frag = fragment ^ (1ULL << atom)
                            # add extended fragments, if not yet present, to the collection
                            if frag not in total_fragments:
                                total_fragments.add(frag)
                                bbsp = self.score_fragment(frag)
                                if bbsp.score < (pars.missingfragmentpenalty + 5):
                                    self.add_fragment(
                                        frag, self.calc_fragment_mass(frag), bbsp.score, bbsp.breaks)
        self.convert_fragments_table()
        return len(self.fragment_info)

    cdef bond_breaks_score_pair score_fragment(self, unsigned long long fragment):
        cdef int b, bondbreaks
        cdef unsigned long long bond
        cdef float score
        cdef bond_breaks_score_pair bbsp
        score = 0
        bondbreaks = 0
        for b in range(self.nbonds):
            bond = self.bonds[b]
            if 0 < (fragment & bond) < bond:
                score += self.bondscore[b]
                bondbreaks += 1
        bbsp.breaks = bondbreaks
        bbsp.score = score
        return bbsp

    def score_fragment_rel2parent(self, unsigned long long fragment, unsigned long long parent):
        cdef int b
        cdef unsigned long long bond
        cdef float score
        score = 0
        for b in range(self.nbonds):
            bond = self.bonds[b]
            if 0 < (fragment & bond) < (bond & parent):
                score += self.bondscore[b]
        return score

    cdef double calc_fragment_mass(self, unsigned long long fragment):
        cdef int atom
        cdef double fragment_mass = 0.0
        for atom in range(self.natoms):
            if fragment & (1ULL << atom):
                fragment_mass += self.atom_masses[atom]
        return fragment_mass

    def add_fragment(self, unsigned long long fragment, double fragmentmass, score, int bondbreaks):
        mass_range = ((self.max_broken_bonds + self.max_water_losses - bondbreaks) * [0] +
                      list(numpy.arange(-bondbreaks + self.ionisation_mode * (1 - self.molcharge),
                                        bondbreaks + self.ionisation_mode * (1 - self.molcharge) + 1) * pars.Hmass + fragmentmass) +
                      (self.max_broken_bonds + self.max_water_losses - bondbreaks) * [0])
        if bondbreaks == 0:
            # make sure that fragmentmass is included
            mass_range[self.max_broken_bonds + self.max_water_losses -
                       self.ionisation_mode] = fragmentmass
        self.fragment_masses += mass_range
        self.fragment_info.append([fragment, score, bondbreaks])

    def convert_fragments_table(self):
        self.fragment_masses_np = numpy.array(self.fragment_masses).reshape(
            len(self.fragment_info), (self.max_broken_bonds + self.max_water_losses) * 2 + 1)

    def calc_avg_score(self):
        self.avg_score = numpy.average(self.scores)

    def get_avg_score(self):
        return self.avg_score

    def find_fragments(self, mass, parent, precision, mz_precision_abs):
        cdef int i
        result = numpy.where(numpy.where(self.fragment_masses_np < max(mass * precision, mass + mz_precision_abs),
                                         self.fragment_masses_np, 0) > min(mass / precision, mass - mz_precision_abs))
        fragment_set = []
        for i in range(len(result[0])):
            fid = result[0][i]
            fragment_set.append(self.fragment_info[fid] +
                                [self.fragment_masses_np[fid][self.max_broken_bonds + self.max_water_losses - self.ionisation_mode * (1 - self.molcharge)]] +
                                [self.ionisation_mode * (1 - self.molcharge) + result[1][i] - self.max_broken_bonds - self.max_water_losses])
        return fragment_set

    def get_fragment_info(self, unsigned long long fragment, deltaH):
        cdef int atom
        mol = Chem.MolFromMolBlock(str(self.mol))
        atomlist = []
        elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'F': 0,
                    'P': 0, 'S': 0, 'Cl': 0, 'Br': 0, 'I': 0}
        for atom in range(self.natoms):
            if ((1ULL << atom) & fragment):
                atomlist.append(atom)
                elements[self.atom_elements[atom]] += 1
                elements['H'] += self.atomHs[atom]
        formula = ''
        for el in ('C', 'H', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'):
            nel = elements[el]
            if nel > 0:
                formula += el
            if nel > 1:
                formula += str(nel)
        atomstring = ','.join(str(a) for a in atomlist)
        return atomstring, atomlist, formula, fragment2inchikey(mol, atomlist)

    def get_natoms(self):
        return self.natoms

    def accepted(self):
        return (self.accept == 1)


def fragment2inchikey(mol, atomlist):
    emol = Chem.EditableMol(mol)
    for atom in reversed(range(mol.GetNumAtoms())):
        if atom not in atomlist:
            emol.RemoveAtom(atom)
    frag = emol.GetMol()
    return Chem.MolToSmiles(frag)
