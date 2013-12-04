#### cython: profile=True
import numpy
cimport numpy
import pars
import ConfigParser, os
config = ConfigParser.ConfigParser()
config.read(['magma_job.ini', os.path.expanduser('~/magma_job.ini')])

if config.get('magma job','chemical_engine')=="rdkit":
    import rdkit_engine as Chem     # Use rdkit_engine
elif config.get('magma job','chemical_engine')=="cdk":
    import cdk_engine               # Use cdk_engine
    Chem=cdk_engine.engine()

ctypedef struct bonded_atom:
    int nbonds
    int[8] atoms

ctypedef struct bond_breaks_score_pair:
    int breaks
    float score

cdef class FragmentEngine(object):

    cdef unsigned long long new_fragment,template_fragment
    cdef int max_broken_bonds,max_water_losses,ionisation_mode
    cdef bonded_atom[64] bonded_atoms
    cdef double[64] atom_masses
    cdef list neutral_loss_atoms
    cdef int nbonds, natoms, accept, skip_fragmentation
    cdef unsigned long long[128] bonds
    cdef float[128] bondscore
    cdef numpy.ndarray fragment_masses_np
    cdef list fragment_masses,fragment_info
    cdef int[64] atomHs
    cdef dict atom_elements
    cdef char* mol
    
    
    def __init__(self,mol,max_broken_bonds,max_water_losses,ionisation_mode,skip_fragmentation):
        cdef unsigned long long bond
        cdef float bondscore
        cdef int x,a1,a2
        
        self.mol=mol
        try:
            mol=Chem.MolFromMolBlock(str(self.mol))
            self.accept=1
            self.natoms=Chem.natoms(mol)  # number of atoms in the molecule
        except:
            self.accept=0
            return
        if self.natoms>64:
            self.accept=0
            return
        self.max_broken_bonds=max_broken_bonds
        self.max_water_losses=max_water_losses
        self.ionisation_mode=ionisation_mode
        self.skip_fragmentation=skip_fragmentation
        self.nbonds=Chem.nbonds(mol)
        self.neutral_loss_atoms=[]
        self.atom_elements={}
        # self.atom_masses=[]
        # self.bonded_atoms=[]           # [[list of atom numbers]]
        # self.bonds=set([])
        # self.bondscore={}
        self.new_fragment=0
        self.template_fragment=0
        self.fragment_masses=((max_broken_bonds+max_water_losses)*2+1)*[0]
        self.fragment_info=[[0,0,0]]
        # self.avg_score=None
        
        for x in range(self.natoms):
            self.bonded_atoms[x].nbonds=0
            self.atom_masses[x]=Chem.GetExtendedAtomMass(mol,x)
            self.atomHs[x]=Chem.GetAtomHs(mol,x)
            self.atom_elements[x]=Chem.GetAtomSymbol(mol,x)
            if Chem.GetAtomSymbol(mol,x) == 'O' and Chem.GetAtomHs(mol,x) == 1 and Chem.GetNBonds(mol,x)==1:
                self.neutral_loss_atoms.append(x)
            if Chem.GetAtomSymbol(mol,x) == 'N' and Chem.GetAtomHs(mol,x) == 2 and Chem.GetNBonds(mol,x)==1:
                self.neutral_loss_atoms.append(x)
        for x in range(self.nbonds):
            a1,a2 = Chem.GetBondAtoms(mol,x)
            self.bonded_atoms[a1].atoms[self.bonded_atoms[a1].nbonds]=a2
            self.bonded_atoms[a1].nbonds+=1
            self.bonded_atoms[a2].atoms[self.bonded_atoms[a2].nbonds]=a1
            self.bonded_atoms[a2].nbonds+=1
            bond = (1ULL<<a1) | (1ULL<<a2)
            bondscore = pars.typew[Chem.GetBondType(mol,x)]*pars.heterow[Chem.GetAtomSymbol(mol,a1) != 'C' or Chem.GetAtomSymbol(mol,a2) != 'C']
            self.bonds[x]=bond
            self.bondscore[x]=bondscore
                
    cdef void extend(self,int atom):
        cdef int a,bonded_a
        cdef unsigned long long atombit
        for a in range(self.bonded_atoms[atom].nbonds):
            bonded_a=self.bonded_atoms[atom].atoms[a]
            atombit=1ULL<<bonded_a
            if atombit & self.template_fragment and not atombit & self.new_fragment:
                self.new_fragment = self.new_fragment | atombit
                self.extend(bonded_a)

    def generate_fragments(self):
        cdef unsigned long long fragment,frag
        cdef int atom,a,bonded_a
        cdef bond_breaks_score_pair bbsp
        cdef set all_fragments,total_fragments,current_fragments,new_fragments
        frag=(1ULL<<self.natoms)-1
        all_fragments=set([frag])
        total_fragments=set([frag])
        current_fragments=set([frag])
        new_fragments=set([frag])
        self.add_fragment(frag,self.calc_fragment_mass(frag),0,0)
        # generate fragments

        if self.skip_fragmentation:
            self.convert_fragments_table()
            return len(self.fragment_info)

        for step in range(self.max_broken_bonds):                    # perform fragmentation for nstep steps
            for fragment in current_fragments:   # loop of all fragments to be fragmented
                for atom in range(self.natoms):       # loop of all atoms
                    if (1ULL<<atom) & fragment:            # in the fragment
                        self.template_fragment=fragment^(1ULL<<atom) # remove the atom
                        list_ext_atoms=set([])
                        extended_fragments=set([])
                        for a in range(self.bonded_atoms[atom].nbonds):              # find all its bonded atoms
                            bonded_a=self.bonded_atoms[atom].atoms[a]
                            if (1ULL<<bonded_a) & self.template_fragment:        # present in the fragment
                                list_ext_atoms.add(bonded_a)
                        if len(list_ext_atoms)==1:                         # in case of one bonded atom, the new fragment
                            extended_fragments.add(self.template_fragment) # is the remainder of the old fragment
                        else:
                            for a in list_ext_atoms:                # otherwise extend all atoms
                                for frag in extended_fragments:     # except when deleted atom is in a ring
                                    if (1ULL<<a) & frag:               # -> previous extended fragment contains
                                        break                       #    already the ext_atom, calculate fragment only once
                                else:
                                    self.new_fragment=1ULL<<a          # extend atom
                                    self.extend(a)
                                    extended_fragments.add(self.new_fragment)
                        for frag in extended_fragments:
                            if frag not in all_fragments:   # add extended fragments if not yet present
                                all_fragments.add(frag)     # to the collection
                                bbsp=self.score_fragment(frag)
                                if bbsp.breaks<=self.max_broken_bonds and bbsp.score < (pars.missingfragmentpenalty+5):
                                    new_fragments.add(frag)
                                    total_fragments.add(frag)
                                    self.add_fragment(frag,self.calc_fragment_mass(frag),bbsp.score,bbsp.breaks)
            current_fragments=new_fragments
            new_fragments=set([])
        for step in range(self.max_water_losses):                    # number of OH losses
            for fi in self.fragment_info:                           # loop of all fragments
                if fi[2]==self.max_broken_bonds+step:               # on which to apply neutral loss rules
                    fragment=fi[0]
                    for atom in self.neutral_loss_atoms:       # loop of all atoms
                        if (1ULL<<atom) & fragment:            # in the fragment
                            frag=fragment^(1ULL<<atom)
                            if frag not in total_fragments:   # add extended fragments if not yet present
                                total_fragments.add(frag)     # to the collection
                                bbsp=self.score_fragment(frag)
                                if bbsp.score < (pars.missingfragmentpenalty+5):
                                    self.add_fragment(frag,self.calc_fragment_mass(frag),bbsp.score,bbsp.breaks)
        self.convert_fragments_table()
        return len(self.fragment_info)

    cdef bond_breaks_score_pair score_fragment(self,unsigned long long fragment):
        cdef int b,bondbreaks
        cdef unsigned long long bond
        cdef float score
        cdef bond_breaks_score_pair bbsp
        score=0
        bondbreaks=0
        for b in range(self.nbonds):
            bond=self.bonds[b]
            if 0 < (fragment & bond) < bond:
                score+=self.bondscore[b]
                bondbreaks+=1
        bbsp.breaks=bondbreaks
        bbsp.score=score
        #if score==0:
        #    print "score=0: ",fragment,bondbreaks
        return bbsp

    def score_fragment_rel2parent(self,unsigned long long fragment,unsigned long long parent):
        cdef int b
        cdef unsigned long long bond
        cdef float score
        score=0
        for b in range(self.nbonds):
            bond=self.bonds[b]
            if 0 < (fragment & bond) < (bond & parent):
                score+=self.bondscore[b]
        return score
    
    cdef double calc_fragment_mass(self, unsigned long long fragment):
        cdef int atom
        cdef double fragment_mass=0.0
        for atom in range(self.natoms):
            if fragment & (1ULL<<atom):
                fragment_mass+=self.atom_masses[atom]
        return fragment_mass

    def add_fragment(self,unsigned long long fragment,double fragmentmass,score,int bondbreaks):
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
        cdef int i
        result=numpy.where(numpy.where(self.fragment_masses_np < max(mass*precision,mass+mz_precision_abs),
                                 self.fragment_masses_np,0) > min(mass/precision,mass-mz_precision_abs))
        fragment_set=[]
        for i in range(len(result[0])):
            fid=result[0][i]
            fragment_set.append(self.fragment_info[fid]+\
                                 [self.fragment_masses_np[fid][self.max_broken_bonds+self.max_water_losses-self.ionisation_mode]]+\
                                 [self.max_broken_bonds+self.max_water_losses-self.ionisation_mode-result[1][i]])
        return fragment_set
    
    def get_fragment_info(self,unsigned long long fragment,deltaH):
        cdef int atom
        mol=Chem.MolFromMolBlock(str(self.mol))
        atomstring=""
        atomlist=[]
        elements={'C':0,'H':0,'N':0,'O':0,'F':0,'P':0,'S':0,'Cl':0,'Br':0,'I':0}
        for atom in range(self.natoms):
            if ((1ULL<<atom) & fragment):
                atomstring+=','+str(atom)
                atomlist.append(atom)
                elements[Chem.GetAtomSymbol(mol,atom)]+=1
                elements['H']+=Chem.GetAtomHs(mol,atom)
        elements['H']-=deltaH
        formula=''
        for el in ('C','H','N','O','F','P','S','Cl','Br','I'):
            nel=elements[el]
            if nel>0:
                formula+=el
            if nel>1:
                formula+=str(nel)
        return atomstring,atomlist,formula,Chem.FragmentToInchiKey(mol,atomlist)
    
    def get_natoms(self):
        return self.natoms
    
    def accepted(self):
        return (self.accept==1)
