from rdkit.Chem import *
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Descriptors

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
bondtype2string = {v:k for k, v in Chem.rdchem.BondType.names.items()}

class newclass(object):
    """
    Additional class
    """
    pass
#    def __init__(self):
#        Chem.__init__(self)
#        print "rdkit_engine available"
#
#    def test(self):
#        print "rdkit_engine available"

def LogP(mol):
    return Chem.Crippen.MolLogP(mol)
def natoms(mol):
    return mol.GetNumAtoms()
def GetExtendedAtomMass(mol,a):
    atom = mol.GetAtomWithIdx(a)
    return mims[atom.GetAtomicNum()]+Hmass*(atom.GetNumImplicitHs()+atom.GetNumExplicitHs())
def GetAtomSymbol(mol,a):
    return mol.GetAtomWithIdx(a).GetSymbol()
def GetAtomHs(mol,a):
    atom = mol.GetAtomWithIdx(a)
    return atom.GetNumImplicitHs()+atom.GetNumExplicitHs()
def nbonds(mol):
    return mol.GetNumBonds()
def GetBondAtoms(mol,b):
    bond=mol.GetBondWithIdx(b)
    return [bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()]
def GetBondType(mol,b):
    bond=mol.GetBondWithIdx(b)
    return bondtype2string[bond.GetBondType()]
def MolToInchiKey(mol):
    # For some reason inchikey fails when bond flag 5 (ISCONJUGATED) is set
    # Make a copy and set all flag 5 values to 0
    return AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))
def FragmentToInchiKey(mol,atomlist):
    emol = Chem.EditableMol(mol)
    for atom in reversed(atomlist):
        emol.RemoveAtom(atom)
    frag = emol.GetMol()
    return Chem.MolToSmiles(frag)
