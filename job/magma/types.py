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
#    def massmatch_rel(self,mim,low,high):
#        for x in range(low,high+1):
#            # if self.mz/me.precision < mim+x*Hmass < self.mz*me.precision:
#            if self.mz/1.000005 < mim+x*Hmass < self.mz*1.000005:
#            # if mim/precision < self.mz-x*Hmass < mim*precision:
#                return x
#        else:
#            return False

class HitType(object):
    def __init__(self,peak,fragment,score,bondbreaks,mass,deltaH):
        self.mz = peak.mz
        self.intensity = peak.intensity
        self.intensity_weight = peak.missing_fragment_score / missingfragmentpenalty
        self.scan = peak.scan
        self.fragment = fragment
        self.score = score
        self.breaks = bondbreaks
        self.mass = mass
        self.deltaH = deltaH
        self.bonds = []
        self.allbonds = 0
        self.besthits=[]
        self.atomstring=''
        self.atomlist=[]
        self.inchikey=""
        #print "childscan",peak.childscan

