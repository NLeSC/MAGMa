import numpy
import pars
import os
from rdkit import Chem

class FragmentEngine(object):
    def __init__(self,mol,max_broken_bonds,max_water_losses,ionisation_mode,skip_fragmentation,molcharge):
        try:
            self.mol=Chem.MolFromMolBlock(str(mol))
            self.accept=True
            self.natoms=self.mol.GetNumAtoms()  # number of atoms in the molecule
        except:
            self.accept=False
            return
        self.max_broken_bonds=max_broken_bonds
        self.max_water_losses=max_water_losses
        self.ionisation_mode=ionisation_mode
        self.skip_fragmentation=skip_fragmentation
        self.molcharge=molcharge
        self.atom_masses=[]
        self.atomHs=[]
        self.neutral_loss_atoms=[]
        self.bonded_atoms=[]           # [[list of atom numbers]]
        self.bonds=set([])
        self.bondscore={}
        self.new_fragment=0
        self.template_fragment=0
        self.fragment_masses=((max_broken_bonds+max_water_losses)*2+1)*[0]
        self.fragment_info=[[0,0,0]]
        self.avg_score=None

        for x in range(self.natoms):
            self.bonded_atoms.append([])
            atom = self.mol.GetAtomWithIdx(x)
            self.atomHs.append(atom.GetNumImplicitHs()+atom.GetNumExplicitHs())
            self.atom_masses.append(pars.mims[atom.GetSymbol()]+pars.Hmass*(self.atomHs[x]))
            if atom.GetSymbol() == 'O' and self.atomHs[x] == 1 and len(atom.GetBonds())==1:
                self.neutral_loss_atoms.append(x)
            if atom.GetSymbol() == 'N' and self.atomHs[x] == 2 and len(atom.GetBonds())==1:
                self.neutral_loss_atoms.append(x)
        for bond in self.mol.GetBonds():
            a1,a2 = bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
            self.bonded_atoms[a1].append(a2)
            self.bonded_atoms[a2].append(a1)
            bondbits = 1<<a1 | 1<<a2
            bondscore = pars.typew[bond.GetBondType()]*pars.heterow[bond.GetBeginAtom().GetSymbol() != 'C' or bond.GetEndAtom().GetSymbol() != 'C']
            self.bonds.add(bondbits)
            self.bondscore[bondbits]=bondscore

    def extend(self,atom):
        for a in self.bonded_atoms[atom]:
            atombit=1<<a
            if atombit & self.template_fragment and not atombit & self.new_fragment:
                self.new_fragment = self.new_fragment | atombit
                self.extend(a)

    def generate_fragments(self):
        frag=(1<<self.natoms)-1
        all_fragments=set([frag])
        total_fragments=set([frag])
        current_fragments=set([frag])
        new_fragments=set([frag])
        self.add_fragment(frag,self.calc_fragment_mass(frag),0,0)
        
        if self.skip_fragmentation:
            self.convert_fragments_table()
            return len(self.fragment_info)
        
        # generate fragments
        for step in range(self.max_broken_bonds):                    # perform fragmentation for nstep steps
            for fragment in current_fragments:   # loop of all fragments to be fragmented
                for atom in range(self.natoms):       # loop of all atoms
                    if (1<<atom) & fragment:            # in the fragment
                        self.template_fragment=fragment^(1<<atom) # remove the atom
                        list_ext_atoms=set([])
                        extended_fragments=set([])
                        for a in self.bonded_atoms[atom]:              # find all its bonded atoms
                            if (1<<a) & self.template_fragment:        # present in the fragment
                                list_ext_atoms.add(a)
                        if len(list_ext_atoms)==1:                         # in case of one bonded atom, the new fragment
                            extended_fragments.add(self.template_fragment) # is the remainder of the old fragment
                        else:
                            for a in list_ext_atoms:                # otherwise extend all atoms
                                for frag in extended_fragments:     # except when deleted atom is in a ring
                                    if (1<<a) & frag:               # -> previous extended fragment contains
                                        break                       #    already the ext_atom, calculate fragment only once
                                else:
                                    self.new_fragment=1<<a          # extend atom
                                    self.extend(a)
                                    extended_fragments.add(self.new_fragment)
                        for frag in extended_fragments:
                            if frag not in all_fragments:   # add extended fragments if not yet present
                                all_fragments.add(frag)     # to the collection
                                bondbreaks,score=self.score_fragment(frag)
                                if bondbreaks<=self.max_broken_bonds and score < (pars.missingfragmentpenalty+5):
                                    new_fragments.add(frag)
                                    total_fragments.add(frag)
                                    self.add_fragment(frag,self.calc_fragment_mass(frag),score,bondbreaks)
            current_fragments=new_fragments
            new_fragments=set([])
        for step in range(self.max_water_losses):                    # number of OH losses
            for fi in self.fragment_info:                           # loop of all fragments
                if fi[2]==self.max_broken_bonds+step:               # on which to apply neutral loss rules
                    fragment=fi[0]
                    for atom in self.neutral_loss_atoms:       # loop of all atoms
                        if (1<<atom) & fragment:            # in the fragment
                            frag=fragment^(1<<atom)
                            if frag not in total_fragments:   # add extended fragments if not yet present
                                total_fragments.add(frag)     # to the collection
                                bondbreaks,score=self.score_fragment(frag)
                                if score < (pars.missingfragmentpenalty+5):
                                    self.add_fragment(frag,self.calc_fragment_mass(frag),score,bondbreaks)
        self.convert_fragments_table()
        return len(self.fragment_info)

    def score_fragment(self,fragment):
        score=0
        bondbreaks=0
        for bond in self.bonds:
            if 0 < (fragment & bond) < bond:
                score+=self.bondscore[bond]
                bondbreaks+=1
        if score==0:
            print "score=0: ",fragment,bondbreaks
        return bondbreaks,score

    def score_fragment_rel2parent(self,fragment,parent):
        score=0
        for bond in self.bonds:
            if 0 < (fragment & bond) < (bond & parent):
                score+=self.bondscore[bond]
        return score
    
    def calc_fragment_mass(self,fragment):
        fragment_mass=0.0
        for atom in range(self.natoms):
            if fragment & (1<<atom):
                fragment_mass+=self.atom_masses[atom]
        return fragment_mass

    def add_fragment(self, fragment, fragmentmass, score, bondbreaks):
        mass_range=((self.max_broken_bonds+self.max_water_losses-bondbreaks)*[0]+\
                  list(numpy.arange(-bondbreaks+self.ionisation_mode*(1-self.molcharge),bondbreaks+self.ionisation_mode*(1-self.molcharge)+1)*pars.Hmass+fragmentmass)+\
                  (self.max_broken_bonds+self.max_water_losses-bondbreaks)*[0])
        if bondbreaks==0:
            mass_range[self.max_broken_bonds+self.max_water_losses-self.ionisation_mode]=fragmentmass # make sure that fragmentmass is included
        self.fragment_masses+=mass_range
        self.fragment_info.append([fragment,score,bondbreaks])
    
    def convert_fragments_table(self):
        self.fragment_masses_np=numpy.array(self.fragment_masses).reshape(len(self.fragment_info),(self.max_broken_bonds+self.max_water_losses)*2+1)

    def calc_avg_score(self):
        # self.avg_score = sum([i[1] for i in self.info])/len(self.info)
        self.avg_score = numpy.average(self.scores)

    def get_avg_score(self):
        return self.avg_score

    def find_fragments(self,mass,parent,precision,mz_precision_abs):
        result=numpy.where(numpy.where(self.fragment_masses_np < max(mass*precision,mass+mz_precision_abs),
                                 self.fragment_masses_np,0) > min(mass/precision,mass-mz_precision_abs))
        fragment_set=[]
        for i in range(len(result[0])):
            fid=result[0][i]
            #if self.max_broken_bonds+self.max_water_losses-self.ionisation_mode-result[1][i] in [0,-self.ionisation_mode]:
            fragment_set.append(self.fragment_info[fid]+\
                                 [self.fragment_masses_np[fid][self.max_broken_bonds+self.max_water_losses-self.ionisation_mode*(1-self.molcharge)]]+\
                                 [self.ionisation_mode*(1-self.molcharge)+result[1][i]-self.max_broken_bonds-self.max_water_losses])
        return fragment_set
    
    def get_fragment_info(self,fragment,deltaH):
        atomstring=""
        atomlist=[]
        elements={'C':0,'H':0,'N':0,'O':0,'F':0,'P':0,'S':0,'Cl':0,'Br':0,'I':0}
        for atom in range(self.natoms):
            if ((1<<atom) & fragment):
                atomstring+=','+str(atom)
                atomlist.append(atom)
                elements[self.mol.GetAtomWithIdx(atom).GetSymbol()]+=1
                elements['H']+=self.atomHs[atom]
        formula=''
        for el in ('C','H','N','O','F','P','S','Cl','Br','I'):
            nel=elements[el]
            if nel>0:
                formula+=el
            if nel>1:
                formula+=str(nel)
        return atomstring,atomlist,formula,fragment2inchikey(self.mol,atomlist)

    def get_natoms(self):
        return self.natoms

    def accepted(self):
        return self.accept

def fragment2inchikey(mol,atomlist):
    emol = Chem.EditableMol(mol)
    for atom in reversed(range(mol.GetNumAtoms())):
        if atom not in atomlist:
            emol.RemoveAtom(atom)
    frag = emol.GetMol()
    return Chem.MolToSmiles(frag)
