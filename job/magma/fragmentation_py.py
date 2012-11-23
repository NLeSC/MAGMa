import rdkit_engine as Chem
import numpy


typew={"AROMATIC":3.0,\
       "DOUBLE":2.0,\
       "TRIPLE":3.0,\
       "SINGLE":1.0}
global missingfragmentpenalty
ringw={False:1,True:1}
heterow={False:2,True:1}
missingfragmentpenalty=10


mims={1:1.0078250321,\
      6:12.0000000,\
      7:14.0030740052,\
      8:15.9949146221,\
      9:18.99840320,\
      15:30.97376151,\
      16:31.97207069,\
      17:34.96885271,\
      35:78.9183376,\
      53:126.904468}

Hmass=mims[1]     # Mass of hydrogen atom


class FragmentEngine(object):
    def __init__(self,structure,max_broken_bonds,max_small_losses):
        self.mol=Chem.MolFromMolBlock(str(structure.mol))
        self.max_broken_bonds=max_broken_bonds
        self.max_small_losses=max_small_losses
        self.natoms=Chem.natoms(self.mol)  # number of atoms in the molecule
        self.atom_masses=[]
        self.neutral_loss_atoms=[]
        self.bonded_atoms=[]           # [[list of atom numbers]]
        self.bonds=set([])
        self.bondscore={}
        self.new_fragment=0
        self.template_fragment=0
        self.fragment_masses=numpy.zeros((max_broken_bonds+max_small_losses)*2+3)
        self.fragments=numpy.array([0])
        self.bondbreaks=numpy.array([0])
        self.scores=numpy.array([0.0])
        self.avg_score=None
        frag=(1<<self.natoms)-1

        for x in range(self.natoms):
            self.bonded_atoms.append([])
            self.atom_masses.append(Chem.GetExtendedAtomMass(self.mol,x))
            if Chem.GetAtomSymbol(self.mol,x) == 'O' and Chem.GetAtomHs(self.mol,x) == 1:
                self.neutral_loss_atoms.append(x)
        for x in range(Chem.nbonds(self.mol)):
            a1,a2 = Chem.GetBondAtoms(self.mol,x)
            self.bonded_atoms[a1].append(a2)
            self.bonded_atoms[a2].append(a1)
            bond = 1<<a1 | 1<<a2
            bondscore = typew[Chem.GetBondType(self.mol,x)]*heterow[Chem.GetAtomSymbol(self.mol,a1) != 'C' or Chem.GetAtomSymbol(self.mol,a2) != 'C']
            self.bonds.add(bond)
            self.bondscore[bond]=bondscore
            
        self.all_fragments=set([frag])
        self.total_fragments=set([frag])
        self.current_fragments=set([frag])
        self.new_fragments=set([frag])
        self.add_fragment(frag,self.calc_fragment_mass(frag),0,0)
    
    def extend(self,atom):
        for a in self.bonded_atoms[atom]:
            atombit=1<<a
            if atombit & self.template_fragment and not atombit & self.new_fragment:
                self.new_fragment = self.new_fragment | atombit
                self.extend(a)

    def generate_fragments(self):
        # generate fragments
        for step in range(self.max_broken_bonds):                    # perform fragmentation for nstep steps
            for fragment in self.current_fragments:   # loop of all fragments to be fragmented
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
                            if frag not in self.all_fragments:   # add extended fragments if not yet present
                                self.all_fragments.add(frag)     # to the collection
                                bondbreaks,score=self.score_fragment(frag)
                                if bondbreaks<=self.max_broken_bonds and score < (missingfragmentpenalty+5):
                                    self.new_fragments.add(frag)
                                    self.total_fragments.add(frag)
                                    self.add_fragment(frag,self.calc_fragment_mass(frag),score,bondbreaks)
            self.current_fragments=self.new_fragments
            self.new_fragments=set([])
        for step in range(self.max_small_losses):                    # number of OH losses
            for fragment in self.current_fragments:   # loop of all fragments on which to apply neutral loss rules
                for atom in self.neutral_loss_atoms:       # loop of all atoms
                    if (1<<atom) & fragment:            # in the fragment
                        frag=fragment^(1<<atom)
                        if frag not in self.total_fragments:   # add extended fragments if not yet present
                            self.total_fragments.add(frag)     # to the collection
                            bondbreaks,score=self.score_fragment(frag)
                            if score < (missingfragmentpenalty+5):
                                self.new_fragments.add(frag)
                                self.add_fragment(frag,self.calc_fragment_mass(frag),score,bondbreaks)
            self.current_fragments=self.new_fragments
            self.new_fragments=set([])

        # calculate masses and scores for fragments
        # first items fragment_masses and fragment_info represent the complete molecule
        self.calc_avg_score()

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

    def add_fragment(self,fragment,fragmentmass,score,bondbreaks):
        self.fragment_masses = numpy.vstack((self.fragment_masses,
                       numpy.hstack((numpy.zeros(self.max_broken_bonds+self.max_small_losses-bondbreaks),
                                  numpy.arange(-bondbreaks-1,bondbreaks+2)*Hmass+fragmentmass,
                                  numpy.zeros(self.max_broken_bonds+self.max_small_losses-bondbreaks)
                                  ))
                       ))
        # self.info.append([fragment,score,bondbreaks])
        self.fragments = numpy.hstack((self.fragments,numpy.array([fragment])))
        self.bondbreaks = numpy.hstack((self.bondbreaks,numpy.array([bondbreaks])))
        self.scores = numpy.hstack((self.scores,numpy.array([score])))

    def calc_avg_score(self):
        # self.avg_score = sum([i[1] for i in self.info])/len(self.info)
        self.avg_score = numpy.average(self.scores)

    def get_avg_score(self):
        return self.avg_score

    def find_fragments(self,mass,parent,precision,mz_precision_abs):
        result=numpy.where(numpy.where(self.fragment_masses < max(mass*precision,mass+mz_precision_abs),
                                 self.fragment_masses,0) > min(mass/precision,mass-mz_precision_abs))
        fragment_set=[]
        for i in range(len(result[0])):
            fid=result[0][i]
            fragment_set.append([self.fragments[fid],
                                 self.scores[fid],
                                 self.bondbreaks[fid],
                                 self.fragment_masses[fid][self.max_broken_bonds+self.max_small_losses+1],
                                 self.max_broken_bonds+self.max_small_losses+1-result[1][i]
                                 ])
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
