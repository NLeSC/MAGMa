#!/usr/bin/env python

import sys,base64,subprocess
import sqlite3,struct
import pkg_resources
import numpy as np
import logging
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Descriptors
from lxml import etree
from sqlalchemy import create_engine,and_,desc
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql import func
from models import Base, Metabolite, Scan, Peak, Fragment, Run

"""
RDkit dependencies:
calculate molecular formula
read smiles
generate 2D conformation for those
generate smiles
"""

typew={Chem.rdchem.BondType.AROMATIC:3.0,\
         Chem.rdchem.BondType.DOUBLE:2.0,\
         Chem.rdchem.BondType.TRIPLE:3.0,\
         Chem.rdchem.BondType.SINGLE:1.0}
# typew={Chem.rdchem.BondType.AROMATIC:6.0,\
#          Chem.rdchem.BondType.DOUBLE:2.0,\
#          Chem.rdchem.BondType.TRIPLE:2.0,\
#          Chem.rdchem.BondType.SINGLE:1.0}
# ringw={False:1.0,True:2.0}
ringw={False:1.0,True:1.0}
heterow={False:1.0,True:0.5}
missingfragmentpenalty=10.0
missingpreservedpenalty=5.0

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
# precision=1.000005 # =5 ppm precision to identify masses
#precision=0.001 # 0.001 m/z precision to identify masses
#Nbonds = 4         # Allowed number of bond breaks TODO move to class
#MSfilter = 2e5    # Intensity cutoff (absolute) to select ions
#MSMSfilter = 0.01 # Intensity cutoff (relative to basePeakIntensity) to select fragment peaks
#MSfilter = 1e6    # Intensity cutoff (absolute) to select ions
#MSMSfilter = 0.5 # Intensity cutoff (relative to basePeakIntensity) to select fragment peaks
#MSfilter = 0.0    # Intensity cutoff (absolute) to select ions
#MSMSfilter = 0.05 # Intensity cutoff (relative to basePeakIntensity) to select fragment peaks
#pp = 0.005          # precision for matching precursor mz
#useFragments=True # Assign fragment data
#ionisation = -1  # positive ionisation mode TODO move to class
#useMSMSonly=True # TODO move to class
#maxMSlevel = 2 # TODO move to class

class MagmaSession(object):
    def __init__(self,db_name,description=""):
        engine = create_engine('sqlite:///'+db_name)
        session = sessionmaker()
        session.configure(bind=engine)
        self.db_session = session()
        Base.metadata.create_all(engine)
        try:
            rundata=self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.description == None:
            rundata.description=description
        self.db_session.add(rundata)
        self.db_session.commit()
    def get_structure_engine(self,
                 metabolism_types="phase1,phase2",
                 n_reaction_steps=2
                 ):
        return StructureEngine(self.db_session,metabolism_types,n_reaction_steps)
    def get_ms_data_engine(self,
                 abs_peak_cutoff=1000,
                 rel_peak_cutoff=0.01,
                 max_ms_level=10
                 ):
        return MsDataEngine(self.db_session,abs_peak_cutoff,rel_peak_cutoff,max_ms_level)
    def get_annotate_engine(self,
                 ionisation_mode=1,
                 skip_fragmentation=False,
                 max_broken_bonds=4,
                 ms_intensity_cutoff=1e6,
                 msms_intensity_cutoff=0.1,
                 mz_precision=0.001,
                 precursor_mz_precision=0.005,
                 use_all_peaks=False
                 ):
        return AnnotateEngine(self.db_session,ionisation_mode,skip_fragmentation,max_broken_bonds,
                 ms_intensity_cutoff,msms_intensity_cutoff,mz_precision,precursor_mz_precision,use_all_peaks)
    def commit(self):
        self.db_session.commit()
    def close(self):
        self.db_session.close()


class StructureEngine(object):
    def __init__(self,db_session,metabolism_types,n_reaction_steps):
        self.db_session = db_session
        try:
            rundata=self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.metabolism_types == None:
            rundata.metabolism_types=metabolism_types
        if rundata.n_reaction_steps == None:
            rundata.n_reaction_steps=n_reaction_steps
        self.db_session.add(rundata)
        self.db_session.commit()
        self.metabolism_types=rundata.metabolism_types.split(',')
        self.n_reaction_steps=rundata.n_reaction_steps

    def calc_mim(self,mol):
        mim=0
        for atom in mol.GetAtoms():
            mim+=(mims[atom.GetAtomicNum()]+Hmass*\
              (atom.GetNumImplicitHs()+atom.GetNumExplicitHs()))
        return mim

    def add_structure(self,mol,name,prob,level,sequence,isquery):
        m=Chem.MolFromMolBlock(mol)
        smiles=Chem.MolToSmiles(m)
        molform=Chem.Descriptors.MolecularFormula(m)
        mim=self.calc_mim(m)
        metab=Metabolite(
            mol=unicode(mol, 'utf-8', 'xmlcharrefreplace'), level=level, probability=prob,
            reactionsequence=sequence, smiles=smiles,
            molformula=molform, isquery=isquery, origin=unicode(name, 'utf-8', 'xmlcharrefreplace'),
            mim=mim,
            logp=Chem.Crippen.MolLogP(m)
            )
        try:
            dupid = self.db_session.query(Metabolite).filter_by(smiles=smiles).one()
            if dupid.probability < prob:
                self.db_session.delete(dupid)
                self.db_session.commit()
                self.db_session.add(metab)
                self.db_session.commit()
                sys.stderr.write('Duplicate structure: '+sequence+' '+smiles+' - old one removed\n')
            else:
                sys.stderr.write('Duplicate structure: '+sequence+' '+smiles+' - kept old one\n')
                metab = dupid
        except NoResultFound:
            self.db_session.add(metab)
            self.db_session.commit()
        finally:
            return metab.metid
            # print 'Added structure:',sequence

    def add_structure_tmp(self,mol,name,prob,level,sequence,isquery):
        m=Chem.MolFromMolBlock(mol)
        smiles=Chem.MolToSmiles(m)
        molform=Chem.Descriptors.MolecularFormula(m)
        mim=self.calc_mim(m)
        dupid = self.db_session.query(Metabolite).filter_by(smiles=smiles).all()
        while len(dupid)>0:
            smiles+='_'
            dupid = self.dbsession.query(Metabolite).filter_by(smiles=smiles).all()
        metab=Metabolite(
            mol=unicode(mol), level=level, probability=prob,
            reactionsequence=sequence, smiles=smiles,
            molformula=molform, isquery=isquery, origin=unicode(name),
            mim=mim
            )
        self.db_session.add(metab)

    def metabolize(self,metid,metabolism,nsteps):
        try:
            parent = self.db_session.query(Metabolite).filter_by(metid=metid).one()
        except:
            print 'Metabolite record ',metid,' does not exist.'
            return
        exec_reactor=pkg_resources.resource_filename( #@UndefinedVariable
                                                      'magma', 'script/reactor')
        exec_reactor+=" -as -f 0.15 -m "+str(nsteps)
        metabolism_files={
            "phase1": pkg_resources.resource_filename( #@UndefinedVariable
                                                       'magma', "data/sygma_rules.phase1.smirks"),
            "phase2": pkg_resources.resource_filename( #@UndefinedVariable
                                                       'magma', "data/sygma_rules.phase2.smirks")
            }
        for m in metabolism.split(','):
            if m in metabolism_files:
                exec_reactor=exec_reactor+" -q "+metabolism_files[m]
        sys.stderr.write(exec_reactor + '\n')

        reactor=subprocess.Popen(exec_reactor, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        reactor.stdin.write(parent.mol+'$$$$\n')
        reactor.stdin.close()

        line=reactor.stdout.readline()
        while line != "":
            name=line
            mol=line
            isquery=0
            while line != 'M  END\n':
                line=reactor.stdout.readline()
                mol+=line
            line=reactor.stdout.readline()
            while line != '$$$$\n' and line != "":
                if line=='> <Probability>\n':
                    prob=float(reactor.stdout.readline())
                elif line=='> <Level>\n':
                    level=int(reactor.stdout.readline())
                elif line=='> <ReactionSequence>\n':
                    sequence=''
                    line=reactor.stdout.readline()
                    while line != '\n':
                        sequence=line+sequence
                        line=reactor.stdout.readline()
                    if sequence=='PARENT\n':
                        isquery=1
                line=reactor.stdout.readline()
            self.add_structure(mol,name,prob,level,sequence,isquery)
            line=reactor.stdout.readline()
        reactor.stdout.close()
        self.db_session.commit()

    def metabolize_all(self,metabolism,nsteps):
        logging.warn('Metabolize all')
        parentids = self.db_session.query(Metabolite.metid).all()
        # print parentids
        for parentid, in parentids:
            self.metabolize(parentid,metabolism,nsteps)


class MsDataEngine(object):
    def __init__(self,db_session,abs_peak_cutoff,rel_peak_cutoff,max_ms_level):
        self.db_session = db_session
        try:
            rundata=self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.abs_peak_cutoff == None:
            rundata.abs_peak_cutoff=abs_peak_cutoff
        if rundata.rel_peak_cutoff == None:
            rundata.rel_peak_cutoff=rel_peak_cutoff
        if rundata.max_ms_level == None:
            rundata.max_ms_level=max_ms_level
        self.db_session.add(rundata)
        self.db_session.commit()
        self.abs_peak_cutoff=rundata.abs_peak_cutoff
        self.rel_peak_cutoff=rundata.rel_peak_cutoff
        self.max_ms_level=rundata.max_ms_level

    def store_mzxml_file(self,mzxml_file):
        logging.warn('Store mzxml file')
        rundata=self.db_session.query(Run).one()
        if rundata.ms_filename == None:
            rundata.ms_filename=mzxml_file
            self.db_session.add(rundata)
            self.ms_filename=mzxml_file
            tree=etree.parse(mzxml_file)
            root=tree.getroot()
            namespace='{'+root.nsmap[None]+'}'
            for mzxmlScan in root.findall(namespace+"msRun/"+namespace+"scan"):
                self.store_mzxml_scan(mzxmlScan,0,namespace)
        else:
            sys.exit('Attempt to read MS data twice')

    def store_mzxml_scan(self,mzxmlScan,precScan,namespace):
        scan=Scan(
            scanid=int(mzxmlScan.attrib['num']),
            mslevel=int(mzxmlScan.attrib['msLevel']),
            rt=float(mzxmlScan.attrib['retentionTime'].strip('PTS'))/60,
            lowmz=float(mzxmlScan.attrib['lowMz']),
            highmz=float(mzxmlScan.attrib['highMz']),
            basepeakmz=float(mzxmlScan.attrib['basePeakMz']),
            basepeakintensity=float(mzxmlScan.attrib['basePeakIntensity']),
            totioncurrent=float(mzxmlScan.attrib['totIonCurrent']),
            precursorscanid=precScan
            )
        for child in mzxmlScan:
            if child.tag == namespace+'precursorMz':
                scan.precursormz=float(child.text)
                scan.precursorintensity=float(child.attrib['precursorIntensity'])
            if child.tag == namespace+'peaks':
                self.store_mzxml_peaks(scan,child.text)
            if child.tag == namespace+'scan' and int(child.attrib['msLevel'])<=self.max_ms_level:
                self.store_mzxml_scan(child,scan.scanid,namespace)
        self.db_session.add(scan)

    def store_mzxml_peaks(self,scan,line):
        dbpeak_cutoff=scan.basepeakintensity*self.rel_peak_cutoff
        if dbpeak_cutoff<self.abs_peak_cutoff:
            dbpeak_cutoff=self.abs_peak_cutoff
        decoded = base64.decodestring(line)
        tmp_size = len(decoded)/4
        unpack_format1 = ">%df" % tmp_size
        unpacked = struct.unpack(unpack_format1,decoded)
        # loop over odd=mz, even=intensity
        for mz, intensity in zip(unpacked[::2], unpacked[1::2]):
            if intensity > dbpeak_cutoff:
                self.db_session.add(Peak(scanid=scan.scanid,mz=mz,intensity=intensity))

    def store_peak_list(self,scanid,rt,precursormz,basepeak,peaklist):
        self.db_session.add(Scan(
            scanid=scanid,
            mslevel=1,
            rt=rt,
            lowmz=precursormz,
            highmz=precursormz,
            basepeakmz=precursormz,
            basepeakintensity=basepeak[1],
            precursorscanid=0
            ))
        self.db_session.add(Peak(scanid=scanid,mz=precursormz,intensity=basepeak[1]))
        self.db_session.add(Scan(
            scanid=scanid+1,
            mslevel=2,
            rt=rt,
            lowmz=peaklist[0][0],
            highmz=peaklist[-1][0],
            basepeakmz=basepeak[0],
            basepeakintensity=basepeak[1],
            precursormz=precursormz,
            precursorintensity=basepeak[1],
            precursorscanid=scanid
            ))
        for peak in peaklist:
            self.db_session.add(Peak(scanid=scanid+1,mz=peak[0],intensity=peak[1]))


class AnnotateEngine(object):
    def __init__(self,db_session,ionisation_mode,skip_fragmentation,max_broken_bonds,
                 ms_intensity_cutoff,msms_intensity_cutoff,mz_precision,precursor_mz_precision,use_all_peaks):
        self.db_session = db_session
        try:
            rundata=self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.ionisation_mode == None:
            rundata.ionisation_mode=ionisation_mode
        if rundata.skip_fragmentation == None:
            rundata.skip_fragmentation=skip_fragmentation
        if rundata.max_broken_bonds == None:
            rundata.max_broken_bonds=max_broken_bonds
        if rundata.ms_intensity_cutoff == None:
            rundata.ms_intensity_cutoff=ms_intensity_cutoff
        if rundata.msms_intensity_cutoff == None:
            rundata.msms_intensity_cutoff=msms_intensity_cutoff
        if rundata.mz_precision == None:
            rundata.mz_precision=mz_precision
        if rundata.precursor_mz_precision == None:
            rundata.precursor_mz_precision=precursor_mz_precision
        if rundata.use_all_peaks == None:
            rundata.use_all_peaks=use_all_peaks
        self.db_session.add(rundata)
        self.db_session.commit()
        self.ionisation_mode=rundata.ionisation_mode
        self.skip_fragmentation=rundata.skip_fragmentation
        self.max_broken_bonds=rundata.max_broken_bonds
        self.ms_intensity_cutoff=rundata.ms_intensity_cutoff
        self.msms_intensity_cutoff=rundata.msms_intensity_cutoff
        self.mz_precision=rundata.mz_precision
        self.precursor_mz_precision=rundata.precursor_mz_precision
        self.use_all_peaks=rundata.use_all_peaks

        self.scans=[]

    def build_spectra(self):
        logging.warn('Build spectra')
        me=self
        class ScanType(object):
            def __init__(self,scan):
                self.peaks=[]
                self.scanid=scan.scanid
                self.mslevel=scan.mslevel
                if scan.mslevel==1:
                    cutoff=me.ms_intensity_cutoff
                else:
                    cutoff=scan.basepeakintensity*me.msms_intensity_cutoff
                peaks=me.db_session.query(Peak).filter(Peak.scanid==scan.scanid).filter(Peak.intensity>=cutoff).all()
                for peak in peaks:
                    self.peaks.append(PeakType(peak.mz,peak.intensity,scan.scanid))
                childscans=me.db_session.query(Scan).filter(Scan.precursorscanid==scan.scanid).all()
                for childscan in childscans:
                    for peak in self.peaks:
                        if -me.precursor_mz_precision < childscan.precursormz-peak.mz < me.precursor_mz_precision:
                            peak.add_child_scan(childscan)
                            break
                    #else:
                    #    if childscan.precursorintensity >= cutoff:
                    #        self.peaks.append(peaktype(childscan.precursormz,childscan.precursorintensity,self.scanid,run))
                    #        self.peaks[-1].add_child_scan(childscan)

        class PeakType(object):
            def __init__(self,mz,intensity,scanid):
                self.mz=mz
                self.intensity=intensity
                self.scan=scanid
                self.childscan=None
                self.missingfragmentscore=missingfragmentpenalty # *intensity**0.5

            def add_child_scan(self,scan):
                self.childscan=ScanType(scan)
                for peak in self.childscan.peaks:
                    self.missingfragmentscore+=peak.missingfragmentscore

            def massmatch_rel(self,mim,low,high):
                for x in range(low,high+1):
                    if self.mz/me.precision < mim+x*Hmass < self.mz*me.precision:
                    # if mim/precision < self.mz-x*Hmass < mim*precision:
                        return x
                else:
                    return False

            def massmatch(self,mim,low,high):
                for x in range(low,high+1):
                    if self.mz-me.mz_precision < mim+x*Hmass < self.mz+me.mz_precision:
                        return x
                else:
                    return False

        for dbscan in self.db_session.query(Scan).filter(Scan.mslevel==1).all():
            self.scans.append(ScanType(dbscan))

    def search_all_structures(self):
        logging.warn('Searching all structures')
        for structure in self.db_session.query(Metabolite).all():
            self.search_structure(structure)

    def search_structure(self,structure):
        FragmentBreaks={}
        FragmentMass={}
        Fragmented=False
        hits=[]                # [hits]
        atombits=[]          # [1,2,4,8,16,....]
        atommass={}          # {atombit:atommass(incl. attached hydrogens)}
        bondbits=[]
        bondscore=[]
        mol=Chem.MolFromMolBlock(str(structure.mol))
        me=self

        class hittype(object):
            def __init__(self,peak,fragment,deltaH):
                self.mass = FragmentMass[fragment]
                self.fragment = fragment
                self.deltaH = deltaH
                self.breaks = FragmentBreaks[fragment]
                self.mz = peak.mz
                self.intensity = peak.intensity
                self.scan = peak.scan
                self.bonds = []
                self.allbonds = 0
                self.score()
                self.besthits=[]
                #print "childscan",peak.childscan
                if peak.childscan!=None:
                    self.score += self.findhits(peak.childscan,fragment)
            def score(self):
                self.score = 0
                for bond in range(len(bondbits)):
                    if 0 < (self.fragment & bondbits[bond]) < bondbits[bond]:
                        self.score += bondscore[bond]
                        self.bonds.append(bond)
                        self.allbonds=self.allbonds|bondbits[bond] # set all bits of broken bond atoms
                # self.score*=self.intensity**0.5
                # self.score*=self.mz
            def findhits(self,childscan,subfrag):
                totalscore=0
                # if childscan.mslevel==5:
                #     print "MS level",childscan.mslevel,len(childscan.peaks)
                for peak in childscan.peaks:
                    # print '-->',peak.mz
                    # print spectrum.scan,peak.mz,peak.intensity
                    besthit=None
                    result=np.where(np.where(fragmims < (peak.mz+me.mz_precision),fragmims,0) > (peak.mz-me.mz_precision))
                    for i in range(len(result[0])):
#                        print "Value: %s"%frags[result[0][i]][result[1][i]]
                        if frags[result[0][i]] & subfrag == frags[result[0][i]]:
                            hit=hittype(peak,frags[result[0][i]],me.max_broken_bonds+1-result[1][i])
                            if besthit==None or besthit.score > hit.score:
                                besthit=hit
                    self.besthits.append(besthit)
                    if besthit==None:
                        totalscore+=peak.missingfragmentscore
                    else:
                        totalscore+=min(besthit.score,peak.missingfragmentscore)
                return totalscore
            def get_fragment(self):
                atomstring=''
                for atom in range(len(atombits)):
                    if (atombits[atom] & self.fragment):
                        atomstring+=','+str(atom)
                return(atomstring[1:])
            def write_fragments(self,metid,parentfragid):
                global fragid
                fragid+=1
                currentFragid=fragid
                me.db_session.add(Fragment(
                    metid=metid,
                    scanid=self.scan,
                    mz=self.mz,
                    mass=self.mass,
                    score=self.score,
                    parentfragid=parentfragid,
                    atoms=self.get_fragment(),
                    deltah=self.deltaH
                    ))
                if len(self.besthits)>0:
                    for hit in self.besthits:
                        if hit != None: # still need to work out how to deal with missed fragments
                            hit.write_fragments(metid,currentFragid)
    
        def grow(fragment):
            NewFragments=[]
            for bond in range(len(bondbits)):
                if 0 < (fragment & bondbits[bond]) < bondbits[bond]:
                    NewFragments.append(fragment|bondbits[bond])
            FragmentBreaks[fragment]=len(NewFragments)
            if len(NewFragments) <= me.max_broken_bonds+3:
                for NewFragment in NewFragments:
                    if NewFragment not in FragmentMass:
                        FragmentMass[NewFragment]=FragmentMass[fragment]+atommass[NewFragment^fragment]
                        grow(NewFragment)
    
        # for peak in self.db_session.query(Peak).filter(Scan.mslevel)==1.filter(Peak.intensity>MSfilter)

        for scan in self.scans:
            for peak in scan.peaks:
                if not ((not self.use_all_peaks) and peak.childscan==None):
                    protonation=self.ionisation_mode-(structure.molformula.find('+')>=0)*1
                    deltaH=peak.massmatch(structure.mim,protonation,protonation)
                    if type(deltaH)==int:
                        if not Fragmented:
                            for x in range(mol.GetNumAtoms()):
                                atombits.append(2**x)
                                an=mol.GetAtomWithIdx(x).GetAtomicNum()
                                atommass[2**x]=mims[an]+\
                                                    (mol.GetAtomWithIdx(x).GetNumImplicitHs()+\
                                                    mol.GetAtomWithIdx(x).GetNumExplicitHs())*Hmass
                            for x in mol.GetBonds():
                                bondbits.append(atombits[x.GetBeginAtomIdx()]|\
                                                     atombits[x.GetEndAtomIdx()])
                                bondscore.append(typew[x.GetBondType()]*ringw[x.IsInRing()]*\
                                                      heterow[x.GetBeginAtom().GetAtomicNum() != 6 or \
                                                                 x.GetEndAtom().GetAtomicNum() != 6])
                            sys.stderr.write('\nMetabolite '+str(structure.metid)+': '+str(structure.origin)+str(structure.reactionsequence)+'\n')
                            sys.stderr.write('Mim: '+str(structure.mim+Hmass))
                            for atom in atombits:
                                FragmentMass[atom]=atommass[atom]
                                grow(atom)
                            fragmims=np.zeros(me.max_broken_bonds*2+3)
                            frags=[0]
                            for fragment in FragmentMass:
                                if FragmentBreaks[fragment]<=me.max_broken_bonds:
#                                    print np.hstack((np.zeros(me.max_broken_bonds-FragmentBreaks[fragment]),
#                                                     np.arange(-FragmentBreaks[fragment],FragmentBreaks[fragment]+1)*Hmass+FragmentMass[fragment],
#                                                     np.zeros(me.max_broken_bonds-FragmentBreaks[fragment])
#                                                     ))
                                    fragmims = np.vstack((fragmims,
                                                       np.hstack((np.zeros(me.max_broken_bonds-FragmentBreaks[fragment]),
                                                                  np.arange(-FragmentBreaks[fragment]-1,FragmentBreaks[fragment]+2)*Hmass+FragmentMass[fragment],
                                                                  np.zeros(me.max_broken_bonds-FragmentBreaks[fragment])
                                                                  ))
                                                       ))
                                    frags.append(fragment)
                                    
                            sys.stderr.write('N fragments created: '+str(len(FragmentMass))+"\n")
                            sys.stderr.write('N fragments kept: '+str(len(fragmims))+"\n")
                            Fragmented=True
                        global fragid
                        fragid=self.db_session.query(func.max(Fragment.fragid)).scalar()
                        if fragid == None:
                            fragid = 0

                        sys.stderr.write('Scan: '+str(scan.scanid)+' - Mz: '+str(peak.mz)+' - ')
                        # storeFragment(metabolite.metid,scan.precursorscanid,scan.precursorpeakmz,2**len(metabolite.atombits)-1,deltaH)
                        hits.append(hittype(peak,2**len(atombits)-1,-protonation))
                        sys.stderr.write('Score: '+str(hits[-1].score)+'\n')
                        hits[-1].write_fragments(structure.metid,0)
        self.db_session.commit()

    def search_structure_nominal(self,structure):
        FragmentBreaks={}
        FragmentMass={}
        Fragments={}
        Fragmented=False
        hits=[]                # [hits]
        atombits=[]          # [1,2,4,8,16,....]
        atommass={}          # {atombit:atommass(incl. attached hydrogens)}
        bondbits=[]
        bondscore=[]
        mol=Chem.MolFromMolBlock(str(structure.mol))
        me=self

        class hittype(object):
            def __init__(self,peak,fragment,deltaH):
                self.mass = FragmentMass[fragment]
                self.fragment = fragment
                self.deltaH = deltaH
                self.breaks = FragmentBreaks[fragment]
                self.mz = peak.mz
                self.intensity = peak.intensity
                self.scan = peak.scan
                self.bonds = []
                self.allbonds = 0
                self.score()
                self.besthits=[]
                #print "childscan",peak.childscan
                if peak.childscan!=None:
                    self.score += self.findhits(peak.childscan,fragment)
            def score(self):
                self.score = 0
                for bond in range(len(bondbits)):
                    if 0 < (self.fragment & bondbits[bond]) < bondbits[bond]:
                        self.score += bondscore[bond]
                        self.bonds.append(bond)
                        self.allbonds=self.allbonds|bondbits[bond] # set all bits of broken bond atoms
                # self.score*=self.intensity**0.5
                # self.score*=self.mz
            def findhits(self,childscan,subfrag):
                totalscore=0
                # if childscan.mslevel==5:
                #     print "MS level",childscan.mslevel,len(childscan.peaks)
                for peak in childscan.peaks:
                    mz=int(round(peak.mz))
                    # print spectrum.scan,peak.mz,peak.intensity
                    besthit=None
                    for x in range(-me.max_broken_bonds-1,me.max_broken_bonds+2):
                        if mz+x in Fragments:
                            for fragment in Fragments[mz+x]:
                                if fragment & subfrag == fragment and (-Fragments[mz+x][fragment]-1 <= x <= Fragments[mz+x][fragment]+1):
                                    hit=hittype(peak,fragment,x)
                                    if besthit==None or besthit.score > hit.score:
                                        besthit=hit
                    self.besthits.append(besthit)
                    if besthit==None:
                        totalscore+=peak.missingfragmentscore
                        # totalscore+=missingfragmentpenalty*peak.mz
                    else:
                        totalscore+=min(besthit.score,peak.missingfragmentscore)
                        # totalscore+=min(besthit.score,missingfragmentpenalty*peak.mz)
                return totalscore
            def get_fragment(self):
                atomstring=''
                for atom in range(len(atombits)):
                    if (atombits[atom] & self.fragment):
                        atomstring+=','+str(atom)
                return(atomstring[1:])
            def write_fragments(self,metid,parentfragid):
                global fragid
                fragid+=1
                currentFragid=fragid
                me.db_session.add(Fragment(
                    metid=metid,
                    scanid=self.scan,
                    mz=self.mz,
                    mass=self.mass,
                    score=self.score,
                    parentfragid=parentfragid,
                    atoms=self.get_fragment(),
                    deltah=self.deltaH
                    ))
                if len(self.besthits)>0:
                    for hit in self.besthits:
                        if hit != None: # still need to work out how to deal with missed fragments
                            hit.write_fragments(metid,currentFragid)
    
        def grow(fragment):
            NewFragments=[]
            for bond in range(len(bondbits)):
                if 0 < (fragment & bondbits[bond]) < bondbits[bond]:
                    NewFragments.append(fragment|bondbits[bond])
            FragmentBreaks[fragment]=len(NewFragments)
            if len(NewFragments) <= me.max_broken_bonds+3:
                for NewFragment in NewFragments:
                    if NewFragment not in FragmentMass:
                        FragmentMass[NewFragment]=FragmentMass[fragment]+atommass[NewFragment^fragment]
                        grow(NewFragment)
    
        # for peak in self.db_session.query(Peak).filter(Scan.mslevel)==1.filter(Peak.intensity>MSfilter)

        for scan in self.scans:
            for peak in scan.peaks:
                if not ((not self.use_all_peaks) and peak.childscan==None):
                    protonation=self.ionisation_mode-(structure.molformula.find('+')>=0)*1
                    deltaH=peak.massmatch(structure.mim,protonation,protonation)
                    if type(deltaH)==int:
                        if not Fragmented:
                            for x in range(mol.GetNumAtoms()):
                                atombits.append(2**x)
                                an=mol.GetAtomWithIdx(x).GetAtomicNum()
                                atommass[2**x]=mims[an]+\
                                                    (mol.GetAtomWithIdx(x).GetNumImplicitHs()+\
                                                    mol.GetAtomWithIdx(x).GetNumExplicitHs())*Hmass
                            for x in mol.GetBonds():
                                bondbits.append(atombits[x.GetBeginAtomIdx()]|\
                                                     atombits[x.GetEndAtomIdx()])
                                bondscore.append(typew[x.GetBondType()]*ringw[x.IsInRing()]*\
                                                      heterow[x.GetBeginAtom().GetAtomicNum() != 6 or \
                                                                 x.GetEndAtom().GetAtomicNum() != 6])
                            sys.stderr.write('\nMetabolite '+str(structure.metid)+': '+str(structure.origin)+str(structure.reactionsequence)+'\n')
                            sys.stderr.write('Mim: '+str(structure.mim+Hmass))
                            for atom in atombits:
                                FragmentMass[atom]=atommass[atom]
                                grow(atom)
                            for fragment in FragmentMass:
                                if FragmentBreaks[fragment]<=me.max_broken_bonds:
                                    if int(round(FragmentMass[fragment])) not in Fragments:
                                        Fragments[int(round(FragmentMass[fragment]))]={}
                                    Fragments[int(round(FragmentMass[fragment]))][fragment]=FragmentBreaks[fragment]
                            sys.stderr.write('N fragments: '+str(len(FragmentMass))+"\n")
                            Fragmented=True
                        global fragid
                        fragid=self.db_session.query(func.max(Fragment.fragid)).scalar()
                        if fragid == None:
                            fragid = 0

                        sys.stderr.write('Scan: '+str(scan.scanid)+' - Mz: '+str(peak.mz)+' - ')
                        # storeFragment(metabolite.metid,scan.precursorscanid,scan.precursorpeakmz,2**len(metabolite.atombits)-1,deltaH)
                        hits.append(hittype(peak,2**len(atombits)-1,protonation))
                        sys.stderr.write('Score: '+str(hits[-1].score)+'\n')
                        hits[-1].write_fragments(structure.metid,0)
        self.db_session.commit()

