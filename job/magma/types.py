import os
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import pars

missingfragmentpenalty=10

class ScanType(object):
    def __init__(self,scanid,mslevel):
        self.peaks=[]
        self.scanid=scanid
        self.mslevel=mslevel
    
class PeakType(object):
    def __init__(self,mz,intensity,scanid,missing_fragment_score):
        self.mz=mz
        self.intensity=intensity
        self.scan=scanid
        self.childscan=None
        self.missing_fragment_score=missing_fragment_score
#    def massmatch_rel(self,mim,low,high):
#        for x in range(low,high+1):
#            # if self.mz/me.precision < mim+x*Hmass < self.mz*me.precision:
#            if self.mz/1.000005 < mim+x*Hmass < self.mz*1.000005:
#            # if mim/precision < self.mz-x*Hmass < mim*precision:
#                return x
#        else:
#            return False

class HitType(object):
    def __init__(self,peak,fragment,score,bondbreaks,mass,ionmass,ion):
        self.mz = peak.mz
        self.intensity = peak.intensity
        self.intensity_weight = peak.missing_fragment_score / missingfragmentpenalty
        self.scan = peak.scan
        self.fragment = fragment
        self.score = score
        self.breaks = bondbreaks
        self.mass = mass
        self.deltaH = ionmass
        self.bonds = []
        self.allbonds = 0
        self.besthits=[]
        self.atomstring=''
        self.atomlist=[]
        self.smiles=""
        self.formula=""
        self.ion=ion
        #print "childscan",peak.childscan

def CalcMIM(mol):
    mim=0.0
    for a in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(a)
        mim+=pars.mims[atom.GetSymbol()]+pars.Hmass*(atom.GetNumImplicitHs()+atom.GetNumExplicitHs())
    return mim

class MoleculeType(object):
    def __init__(self,molblock,name,refscore,predicted=0,mim=None,natoms=None,inchikey14=None,molform=None,reference=None,logp=None):
        if inchikey14==None or mim==None or molform==None or logp==None or natoms==None:
            mol=Chem.MolFromMolBlock(molblock)
            inchikey14=Chem.AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))[:14]
            mim=CalcMIM(mol)
            molform=Chem.rdMolDescriptors.CalcMolFormula(mol)
            natoms=mol.GetNumHeavyAtoms()
            logp = Chem.Crippen.MolLogP(mol)

        self.molblock = molblock #: molfile as string
        self.refscore = refscore
        self.inchikey14 = inchikey14 #: Smile string
        self.formula = molform #: Molecular formula
        self.predicted = predicted #: Whether metabolite was given as query or is a result a of reaction
        self.name = name #: Name of molecule
        # self.nhits = Column(Integer)
        self.mim = mim #: Monoisotopic mass
        self.natoms = natoms #: Number of non-hydrogen atoms
        self.logp = logp #: Calculated logP
        self.reference = reference

