import os

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
