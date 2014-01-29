import os
import rdkit_engine as Chem

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
        self.inchikey=""
        self.formula=""
        self.ion=ion
        #print "childscan",peak.childscan

class MoleculeType(object):
    def __init__(self,molblock,name,prob,level,isquery=1,mim=None,natoms=None,inchikey=None,molform=None,reference=None,logp=None):
        if inchikey==None or mim==None or molform==None or logp==None or natoms==None:
            mol=Chem.MolFromMolBlock(molblock)
            inchikey=Chem.MolToInchiKey(mol)[:14]
            # inchikey=Chem.MolToSmiles(mol)
            mim,molform=Chem.GetFormulaProps(mol)
            natoms=mol.GetNumHeavyAtoms()
            logp = Chem.LogP(mol)

        self.molblock = molblock #: molfile as string
        self.level = level
        self.probability = prob
        self.inchikey = inchikey #: Smile string
        self.molformula = molform #: Molecular formula
        self.isquery = isquery #: Whether metabolite was given as query or is a result a of reaction
        self.name = name #: Name of molecule
        # self.nhits = Column(Integer)
        self.mim = mim #: Monoisotopic mass
        self.natoms = natoms #: Number of non-hydrogen atoms
        self.logp = logp #: Calculated logP
        self.reference = reference

