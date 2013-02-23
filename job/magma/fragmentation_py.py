import rdkit_engine as Chem
# import cdk_engine               # Use cdk_engine
# Chem=cdk_engine.engine()
import numpy
import pars


class FragmentEngine(object):
    def __init__(self,mol,max_broken_bonds,max_water_losses,ionisation_mode):
        self.mol=Chem.MolFromMolBlock(str(mol))
        self.max_broken_bonds=max_broken_bonds
        self.max_water_losses=max_water_losses
        self.ionisation_mode=ionisation_mode
        self.natoms=Chem.natoms(self.mol)  # number of atoms in the molecule
        self.atom_masses=[]
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
            self.atom_masses.append(Chem.GetExtendedAtomMass(self.mol,x))
            if Chem.GetAtomSymbol(self.mol,x) == 'O' and Chem.GetAtomHs(self.mol,x) == 1 and Chem.GetNBonds(self.mol,x)==1:
                self.neutral_loss_atoms.append(x)
        for x in range(Chem.nbonds(self.mol)):
            a1,a2 = Chem.GetBondAtoms(self.mol,x)
            self.bonded_atoms[a1].append(a2)
            self.bonded_atoms[a2].append(a1)
            bond = 1<<a1 | 1<<a2
            bondscore = pars.typew[Chem.GetBondType(self.mol,x)]*pars.heterow[Chem.GetAtomSymbol(self.mol,a1) != 'C' or Chem.GetAtomSymbol(self.mol,a2) != 'C']
            self.bonds.add(bond)
            self.bondscore[bond]=bondscore
    
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
        self.fragment_masses+=((self.max_broken_bonds+self.max_water_losses-bondbreaks)*[0]+\
                                  list(numpy.arange(-bondbreaks+self.ionisation_mode,bondbreaks+self.ionisation_mode+1)*pars.Hmass+fragmentmass)+\
                                  (self.max_broken_bonds+self.max_water_losses-bondbreaks)*[0])
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
            fragment_set.append(self.fragment_info[fid]+\
                                 [self.fragment_masses_np[fid][self.max_broken_bonds+self.max_water_losses-self.ionisation_mode]]+\
                                 [self.max_broken_bonds+self.max_water_losses-self.ionisation_mode-result[1][i]])
        return fragment_set
    
    def get_fragment_info(self,fragment):
        atomstring=""
        atomlist=[]
        for atom in range(self.natoms):
            if ((1<<atom) & fragment):
                atomstring+=','+str(atom)
                atomlist.append(atom)
        return atomstring,atomlist # ,Chem.FragmentToInchiKey(self.mol,atomlist)

    def get_natoms(self):
        return self.natoms

    def accepted(self):
        return True
