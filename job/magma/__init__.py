#!/usr/bin/env python

import sys,base64,subprocess,StringIO,time,re
import sqlite3,struct,zlib,gzip,copy
import pkg_resources
import numpy
import logging
from lxml import etree
from sqlalchemy import create_engine,and_,desc
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql import func
from models import Base, Metabolite, Scan, Peak, Fragment, Run
import pp
import cPickle as pickle
import types
import pars
import rdkit_engine as Chem     # Use rdkit_engine
# import cdk_engine               # Use cdk_engine
# Chem=cdk_engine.engine()

"""
RDkit dependencies:
calculate molecular formula
read smiles
generate 2D conformation for those
generate smiles
"""

max_small_losses=1


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
                 mz_precision=5.0,
                 mz_precision_abs=0.001,
                 precursor_mz_precision=0.005,
                 use_all_peaks=False
                 ):
        return AnnotateEngine(self.db_session,ionisation_mode,skip_fragmentation,max_broken_bonds,
                 ms_intensity_cutoff,msms_intensity_cutoff,mz_precision,mz_precision_abs,precursor_mz_precision,use_all_peaks)
    def get_data_analysis_engine(self):
        return DataAnalysisEngine(self.db_session)
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

    def add_structure(self,molblock,name,prob,level,sequence,isquery,mass_filter=9999,mim=None,inchikey=None,molform=None,reference=None,logP=None,check_duplicates=False):
        if inchikey==None or mim==None or molform==None or logP==None:
            mol=Chem.MolFromMolBlock(molblock)
        if inchikey == None:
            inchikey=Chem.MolToInchiKey(mol)[:14]
            # inchikey=Chem.MolToSmiles(mol)
        if mim == None or molform == None:
            mim,molform=Chem.GetFormulaProps(mol)
        if mim > mass_filter:
            return
        if logP == None:
            logP = Chem.LogP(mol)
        # mim=self.calc_mim(mol)
        metab=Metabolite(
            mol=unicode(molblock, 'utf-8', 'xmlcharrefreplace'),
            level=level,
            probability=prob,
            reactionsequence=sequence,
            smiles=inchikey,
            molformula=molform,
            isquery=isquery,
            origin=unicode(name, 'utf-8', 'xmlcharrefreplace'),
            nhits=0,
            mim=mim,
            reference=reference,
            logp=logP
            )
        if check_duplicates and len(self.db_session.query(Metabolite).filter_by(smiles=inchikey).all())>0:
#            if dupid.probability < prob:
#                self.db_session.delete(dupid)
#                metab.metid=dupid.metid
#                self.db_session.add(metab)
#                sys.stderr.write('Duplicate structure: '+sequence+' '+inchikey+' - old one removed\n')
#                # TODO remove any fragments related to this structure as well
#            else:
            sys.stderr.write('Duplicate structure: '+sequence+' '+inchikey+' - kept old one\n')
            return
        else:
            self.db_session.add(metab)
            #sys.stderr.write('Added: '+name+'\n')
            self.db_session.flush()
            return metab.metid

    def add_structure_tmp(self,mol,name,prob,level,sequence,isquery):
        m=Chem.MolFromMolBlock(mol)
        smiles=Chem.MolToSmiles(m)
        molform=Chem.Descriptors.MolecularFormula(m)
        mim=self.calc_mim(m)
        dupid = self.db_session.query(Metabolite).filter_by(smiles=smiles).all()
        while len(dupid)>0:
            smiles+='_'
            dupid = self.db_session.query(Metabolite).filter_by(smiles=smiles).all()
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
                                                       'magma', "data/sygma_rules.phase2.smirks"),
            "gut": pkg_resources.resource_filename( #@UndefinedVariable
                                                       'magma', "data/gut.smirks")
            }
        for m in metabolism.split(','):
            if m in metabolism_files:
                exec_reactor=exec_reactor+" -q "+metabolism_files[m]
        sys.stderr.write(exec_reactor + '\n')

        reactor=subprocess.Popen(exec_reactor, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        reactor.stdin.write(parent.mol+'$$$$\n')
        reactor.stdin.close()

        metids=[]
        line=reactor.stdout.readline()
        while line != "":
            name=line[:-1]
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
                        sequence=''
                    reactionsequence=parent.reactionsequence+sequence
                line=reactor.stdout.readline()
            metids.append(self.add_structure(mol,name,prob,level,reactionsequence,isquery))
            line=reactor.stdout.readline()
        reactor.stdout.close()
        self.db_session.commit()
        return metids

    def metabolize_all(self,metabolism,nsteps):
        logging.warn('Metabolize all')
        parentids = self.db_session.query(Metabolite.metid).all()
        # print parentids
        metids=[]
        for parentid, in parentids:
            metids.extend(self.metabolize(parentid,metabolism,nsteps))
        return set(metids)
    
    def retrieve_structures(self,mass):
        dbfilename = '/home/ridderl/chebi/ChEBI_complete_3star.sqlite'
        conn = sqlite3.connect(dbfilename)
        c = conn.cursor()
        result = c.execute('SELECT * FROM molecules WHERE mim BETWEEN ? AND ?' , (mass-0.01,mass+0.01))
        for (id,mim,molblock,smiles,chebi_name) in result:
            self.add_structure(zlib.decompress(molblock),str(chebi_name),1.0,1,"",1)

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
            # mzxmlScan = root.findall(namespace+"msRun/"+namespace+"scan")[0]
                self.store_mzxml_scan(mzxmlScan,0,namespace)
        else:
            sys.exit('Attempt to read MS data twice')
        self.db_session.commit()

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
                if child.attrib.get('precursorScanNum') != None:
                    scan.precursorscanid=float(child.attrib['precursorScanNum'])
            if child.tag == namespace+'peaks':
                self.store_mzxml_peaks(scan,child.text)
            if child.tag == namespace+'scan' and int(child.attrib['msLevel'])<=self.max_ms_level:
#                if scan.mslevel == 1:
#                    cutoff=0.0
#                else:
#                    cutoff=max(scan.basepeakintensity*self.rel_peak_cutoff,self.abs_peak_cutoff)
#                if float(child.attrib['basePeakIntensity'])>cutoff:
                self.store_mzxml_scan(child,scan.scanid,namespace)
        self.db_session.add(scan)
        self.db_session.flush()

    def store_mzxml_peaks(self,scan,line):
#        dbpeak_cutoff=scan.basepeakintensity*self.rel_peak_cutoff
#        if dbpeak_cutoff<self.abs_peak_cutoff:
#            dbpeak_cutoff=self.abs_peak_cutoff
        decoded = base64.decodestring(line)
        tmp_size = len(decoded)/4
        unpack_format1 = ">%df" % tmp_size
        unpacked = struct.unpack(unpack_format1,decoded)
        for mz, intensity in zip(unpacked[::2], unpacked[1::2]):
            if intensity > self.abs_peak_cutoff:
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
        self.db_session.commit()

    def store_manual_tree(self,manual_tree):
        tree_string=open(manual_tree).read()
        tree_list=re.split('([\,\(\)])',tree_string.replace(" ","").replace("\n",""))
        print tree_list
        self.global_scanid = 1
        self.store_manual_subtree(tree_list,0,0,0,1)
    
    def store_manual_subtree(self,tree_list,precursor_scanid,precursor_mz,precursor_intensity,mslevel):
        lowmz=None
        highmz=None
        basepeakmz=None
        basepeakintensity=None
        scanid=self.global_scanid
        npeaks=0
        while len(tree_list)>0 and tree_list[0]!=')':
            #print tree_list[0]
            if tree_list[0].find(':')>=0:
                mz,intensity=tree_list.pop(0).split(':')
                self.db_session.add(Peak(scanid=scanid,mz=mz,intensity=intensity))
                npeaks+=1
                if lowmz==None or mz<lowmz:
                    lowmz=mz
                if highmz==None or mz>highmz:
                    highmz=mz
                if basepeakintensity==None or intensity>basepeakintensity:
                    basepeakmz=mz
                    basepeakintensity=intensity
            if tree_list[0]=='(':
                tree_list.pop(0)
                self.global_scanid+=1
                self.store_manual_subtree(tree_list,scanid,mz,intensity,mslevel+1)
            if tree_list[0]==',' or tree_list[0]=='':
                tree_list.pop(0)
        if npeaks>0:
            self.db_session.add(Scan(
                scanid=scanid,
                mslevel=mslevel,
                lowmz=lowmz,
                highmz=highmz,
                basepeakmz=basepeakmz,
                basepeakintensity=basepeakintensity,
                precursorscanid=precursor_scanid,
                precursormz=precursor_mz,
                precursorintensity=precursor_intensity
                ))
        self.db_session.commit()
        if len(tree_list)>0:
            tree_list.pop(0)
            

class AnnotateEngine(object):
    def __init__(self,db_session,ionisation_mode,skip_fragmentation,max_broken_bonds,
                 ms_intensity_cutoff,msms_intensity_cutoff,mz_precision,mz_precision_abs,
                 precursor_mz_precision,use_all_peaks):
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
        if rundata.mz_precision_abs == None:
            rundata.mz_precision_abs=mz_precision_abs
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
        self.precision=1+rundata.mz_precision/1e6
        self.mz_precision_abs=rundata.mz_precision_abs
        self.precursor_mz_precision=rundata.precursor_mz_precision
        self.use_all_peaks=rundata.use_all_peaks

        self.scans=[]

    def build_spectrum(self,dbscan):
        scan=types.ScanType(dbscan.scanid,dbscan.mslevel)
        if scan.mslevel==1:
            cutoff=self.ms_intensity_cutoff
        else:
            cutoff=dbscan.basepeakintensity*self.msms_intensity_cutoff
        dbpeaks=self.db_session.query(Peak).filter(Peak.scanid==scan.scanid).filter(Peak.intensity>=cutoff).all()
        for dbpeak in dbpeaks:
            scan.peaks.append(types.PeakType(dbpeak.mz,dbpeak.intensity,scan.scanid,pars.missingfragmentpenalty*(dbpeak.intensity**0.5)))
        dbchildscans=self.db_session.query(Scan).filter(Scan.precursorscanid==scan.scanid).all()
        for dbchildscan in dbchildscans:
            # find the highest peak that qualifies as precursor peak for the child spectrum
            prec_intensity=0.0
            for peak in scan.peaks:
                if peak.intensity>prec_intensity and -self.precursor_mz_precision < dbchildscan.precursormz-peak.mz < self.precursor_mz_precision:
                    prec_peak=peak
                    prec_intensity=peak.intensity
            # if present, process the childscan as its child spectrum
            if prec_intensity>0.0:
                prec_peak.childscan=self.build_spectrum(dbchildscan)
                for childpeak in prec_peak.childscan.peaks:
                    prec_peak.missing_fragment_score+=childpeak.missing_fragment_score
#            else:
#                if dbchildscan.precursorintensity >= cutoff:
#                    scan.peaks.append(types.PeakType(dbchildscan.precursormz,dbchildscan.precursorintensity,scan.scanid,missingfragmentpenalty*(dbchildscan.precursorintensity**0.5)))
#                    scan.peaks[-1].childscan=self.build_spectrum(dbchildscan)
#                    for childpeak in scan.peaks[-1].childscan.peaks:
#                        scan.peaks[-1].missing_fragment_score+=childpeak.missing_fragment_score
        return scan
    
    def build_spectra(self,scans='all'):
        logging.warn('Build spectra')
        if scans=='all':
            for dbscan in self.db_session.query(Scan).filter(Scan.mslevel==1).all():
                self.scans.append(self.build_spectrum(dbscan))
        else:
            logging.warn('for scans: '+str(scans))
            for dbscan in self.db_session.query(Scan).filter(Scan.mslevel==1).filter(Scan.scanid.in_(scans)).all():
                self.scans.append(self.build_spectrum(dbscan))
        self.indexed_peaks={}   # dictionary: sets of peaks for each integer m/z value
        for scan in self.scans:
            for peak in scan.peaks:
                if not ((not self.use_all_peaks) and peak.childscan==None):
                    int_mass=int(round(peak.mz))
                    if int_mass not in self.indexed_peaks:
                        self.indexed_peaks[int_mass]=set([])
                    self.indexed_peaks[int_mass].add(peak)
                    print self.write_peak(peak)

    def write_tree(self,scanid):
        for scan in self.scans:
            if scan.scanid==scanid:
                for peak in scan.peaks:
                    if not ((not self.use_all_peaks) and peak.childscan==None):
                        self.write_peak(peak)

    def write_peak(self,peak):
        peak_string="%.6f: %i" % (peak.mz,peak.intensity)
        if peak.childscan!=None:
            peak_string+=' ('
            n=0
            for childpeak in peak.childscan.peaks:
                if n>0:
                    peak_string+=", "
                peak_string+=self.write_peak(childpeak)
                n+=1
            peak_string+=')'
        return peak_string

    def get_chebi_candidates(self):
        dbfilename = '/home/ridderl/chebi/ChEBI_complete_3star.sqlite'
        conn = sqlite3.connect(dbfilename)
        c = conn.cursor()
        db_candidates={}
        # First a dictionary is created with candidate molecules based on all level1 peaks
        # In this way duplicates originating from repeated detection of the same component
        # are removed before attempting to add the candidates to the database
        for scan in self.scans:
            for peak in scan.peaks:
                mass=peak.mz-self.ionisation_mode*pars.Hmass
                if not ((not self.use_all_peaks) and peak.childscan==None):
                    result = c.execute('SELECT * FROM molecules WHERE mim BETWEEN ? AND ?' , (mass/self.precision,mass*self.precision))
                    for (id,mim,molblock,chebi_name) in result:
                        db_candidates[id]=[molblock,chebi_name]
                    print str(mass)+' --> '+str(len(db_candidates))+' candidates'
        return db_candidates

    def get_pubchem_candidates(self,fast,dbfilename='',min_refscore='',max_mz=''):
        where=''
        if dbfilename=='':
            dbfilename='/media/PubChem/Pubchem_MAGMa.db'
        if min_refscore!='':
            where += ' AND refscore >= '+min_refscore
        if fast:
            where += ' AND natoms <= 64' # fast calculations are based on 64 long int, larger molecules are not allowed
        if max_mz=='':
            max_mz='9999'
        mmz=float(max_mz)
        
        conn = sqlite3.connect(dbfilename)
        conn.text_factory=str
        c = conn.cursor()
        struct_engine = StructureEngine(self.db_session,"",0)

        # build sorted list of query masses
        mzs=[]
        for scan in self.scans:
            for peak in scan.peaks:
                if not ((not self.use_all_peaks) and peak.childscan==None) and peak.mz <= mmz:
                    mzs.append(peak.mz)
        mzs.sort()
        # build non-overlapping set of queries around these masses
        queries=[[0,0]]
        for mz in mzs:
            ql=int(1e6*(mz/self.precision-self.ionisation_mode*(pars.Hmass-pars.elmass)))
            qh=int(1e6*(mz*self.precision-self.ionisation_mode*(pars.Hmass-pars.elmass)))
            if queries[-1][0] <= ql <= queries[-1][1]:
                queries[-1][1]=qh
            else:
                queries.append([ql,qh])

        # in case of an empty database, no check for existing duplicates needed
        check_duplicates = (self.db_session.query(Metabolite.metid).count() > 0)
        print 'check: ',check_duplicates,fast

        # All candidates are stored in dbsession, resulting metids are returned
        metids=set([])
        for ql,qh in queries:
            result = c.execute('SELECT * FROM molecules WHERE mim BETWEEN ? AND ? %s' % where, (ql,qh))
            for (cid,mim,natoms,molblock,inchikey,molform,name,refscore,logp) in result:
                metid=struct_engine.add_structure(molblock=zlib.decompress(molblock),
                               name=name+' ('+str(cid)+')',
                               mim=float(mim/1e6),
                               molform=molform,
                               inchikey=inchikey,
                               prob=refscore,
                               level=1,
                               sequence="",
                               isquery=1,
                               reference='<a href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&cmd=Link&LinkName=pccompound_pccompound_sameisotopic_pulldown&from_uid='+\
                                         str(cid)+'">'+str(cid)+' (PubChem)</a>',
                               logP=float(logp)/10.0,
                               check_duplicates=check_duplicates
                               )
                metids.add(metid)
            print str(ql)+','+str(qh)+' --> '+str(len(metids))+' candidates'
        self.db_session.commit()
        return metids

    def search_structures(self,metids=None,ncpus=1,fast=False):
        logging.warn('Searching structures')
        if fast:
            fragmentation_module='magma.fragmentation_cy'
        else:
            fragmentation_module='magma.fragmentation_py'
        if metids==None:
            structures = self.db_session.query(Metabolite).all()
        else:  # split metids in chunks of 500 to avoid error in db_session.query
            structures = []
            while len(metids)>0:
                ids=set([])
                while len(ids)<500 and len(metids)>0:
                    ids.add(metids.pop())
                structures = structures+self.db_session.query(Metabolite).filter(Metabolite.metid.in_(ids)).all()
        global fragid
        fragid=self.db_session.query(func.max(Fragment.fragid)).scalar()
        if fragid == None:
            fragid = 0
        ppservers = ()
        logging.warn('calculating on '+str(ncpus)+' cpus !!!')
        job_server = pp.Server(ncpus, ppservers=ppservers)
        jobs=[]
        for structure in structures:
            # collect all peaks with masses within 3 Da range
            int_mass=int(round(structure.mim+self.ionisation_mode*pars.Hmass))
            try:
                peaks=self.indexed_peaks[int_mass]
            except:
                peaks=set([])
            try:
                peaks=peaks.union(self.indexed_peaks[int_mass-1])
            except:
                pass
            try:
                peaks=peaks.union(self.indexed_peaks[int_mass+1])
            except:
                pass
            if len(peaks)>0:
                jobs.append((structure,
                   job_server.submit(search_structure,(structure.mol,
                              structure.mim,
                              structure.molformula,
                              peaks,
                              self.max_broken_bonds,
                              max_small_losses,
                              self.precision,
                              self.mz_precision_abs,
                              self.use_all_peaks,
                              self.ionisation_mode,
                              fast
                              ),(),(
                              "magma.types",
                              "magma.pars",
                              fragmentation_module
                              )
                           )))
        for structure,job in jobs:
            raw_result=job(raw_result=True)
            hits,sout = pickle.loads(raw_result)
            #print sout
            sys.stderr.write('Metabolite '+str(structure.metid)+': '+str(structure.origin)+'\n')
            structure.nhits=len(hits)
            self.db_session.add(structure)
            for hit in hits:
                sys.stderr.write('Scan: '+str(hit.scan)+' - Mz: '+str(hit.mz)+' - ')
                # storeFragment(metabolite.metid,scan.precursorscanid,scan.precursorpeakmz,2**len(metabolite.atombits)-1,deltaH)
                sys.stderr.write('Score: '+str(hit.score)+'\n')
                # outfile.write("\t"+str(hit.score/fragment_store.get_avg_score()))
                self.store_hit(hit,structure.metid,0)
            self.db_session.flush()
        self.db_session.commit()

    def store_hit(self,hit,metid,parentfragid):
        global fragid
        fragid+=1
        currentFragid=fragid
        score=hit.score
        deltappm=None
        if score != None:
            score=score/hit.intensity_weight
            deltappm=(hit.mz+hit.deltaH*pars.Hmass-hit.mass+self.ionisation_mode*pars.elmass)/hit.mz*1e6
        # print atomlist, Chem.FragmentSmiles(mol,atomlist)
        self.db_session.add(Fragment(
            metid=metid,
            scanid=hit.scan,
            mz=hit.mz,
            mass=hit.mass,
            score=score,
            parentfragid=parentfragid,
            atoms=hit.atomstring,
            inchikey=hit.inchikey,
            deltah=hit.deltaH,
            deltappm=deltappm
            ))
        if len(hit.besthits)>0:
            for childhit in hit.besthits:
                if childhit != None: # still need to work out how to deal with missed fragments
                    self.store_hit(childhit,metid,currentFragid)


class DataAnalysisEngine(object):
    def __init__(self,db_session):
        self.db_session = db_session

    def get_scores(self,scanid):
        return self.db_session.query(Fragment.score,Metabolite.reactionsequence,Metabolite.molformula).\
            join((Metabolite,and_(Fragment.metid==Metabolite.metid))).\
            filter(Fragment.parentfragid==0).\
            filter(Fragment.scanid==scanid).\
            all()
            
    def get_num_peaks(self,scanid):
        return self.db_session.query(Peak).filter(Peak.scanid==scanid).count()
            
    def export_assigned_molecules(self,name):
        for metabolite,peak in self.db_session.query(Metabolite,Peak).filter(Metabolite.metid==Peak.assigned_metid):
            print metabolite.origin.splitlines()[0]
            print metabolite.mol[metabolite.mol.find("\n")+1:-1]
            print "> <ScanID>\n"+str(peak.scanid)+"\n"
            print "> <mz>\n"+str(peak.mz)+"\n"
            print "> <intensity>\n"+str(peak.intensity)+"\n"
            print "> <rt>\n"+str(self.db_session.query(Scan.rt).filter(Scan.scanid==peak.scanid).all()[0][0])+"\n"
            print "> <molecular formula>\n"+metabolite.molformula+"\n"
            print "$$$$"

    def write_SDF(self,molecules=None,columns=None,sortcolumn=None,descend=False):
        if molecules==None:
            if descend:
                molecules=self.db_session.query(Metabolite).order_by(desc(sortcolumn)).all()
            else:
                molecules=self.db_session.query(Metabolite).order_by(sortcolumn).all()
        for molecule in molecules:
            print molecule.mol[:-1]
            if columns==None:
                columns=dir(molecule)
            for column in columns:
                if column[:1] != '_' and column != 'mol' and column != 'metadata':
                    print '> <'+column+'>\n'+str(molecule.__getattribute__(column))+'\n'
            print '$$$$'

def search_structure(mol,mim,molformula,peaks,max_broken_bonds,max_small_losses,precision,mz_precision_abs,use_all_peaks,ionisation_mode,fast):
    # Chem=magma.cdk_engine.engine()      # Use cdk_engine
    Chem=magma.rdkit_engine             # Use rdkit_engine
    pars=magma.pars
    if fast:
        Fragmentation=magma.fragmentation_cy
    else:
        Fragmentation=magma.fragmentation_py

    def massmatch(peak,mim,low,high):
        for x in range(low,high+1):
            #if self.mz-me.mz_precision < mim+x*Hmass < self.mz+me.mz_precision:
            if peak.mz/precision <= mim+x*pars.Hmass-ionisation_mode*pars.elmass <= peak.mz*precision:
                return x
        else:
            return False

    #def findhit(self,childscan,parent):
    def gethit (peak,fragment,score,bondbreaks,mass,deltaH):
        try:
            hit=types.HitType(peak,fragment,score,bondbreaks,mass,deltaH)
        except:
            hit=magma.types.HitType(peak,fragment,score,bondbreaks,mass,deltaH)
        if fragment>0 and peak.childscan!=None and len(peak.childscan.peaks) > 0: # fragment=0 means it is a missing fragment
            n_child_peaks=len(peak.childscan.peaks)
            total_score=0.0
            total_count=0.0
            for childpeak in peak.childscan.peaks:
                besthit=gethit(childpeak,0,None,0,0,0)
                mz_neutral=childpeak.mz+ionisation_mode*pars.elmass # m/z value of the neutral form of the fragment (mass of electron added/removed)
                for childfrag,childscore,childbbreaks,childmass,childH in fragment_engine.find_fragments(mz_neutral,fragment,precision,mz_precision_abs):
                    if childfrag & fragment == childfrag:
                        childhit=gethit(childpeak,childfrag,childscore*(childpeak.intensity**0.5),childbbreaks,childmass,childH)
                        if besthit.score==None or besthit.score > childhit.score or \
                               (besthit.score == childhit.score and abs(besthit.deltaH) > abs(childhit.deltaH)) or \
                               fragment_engine.score_fragment_rel2parent(besthit.fragment,fragment) > fragment_engine.score_fragment_rel2parent(childhit.fragment,fragment):
                            besthit=childhit
                hit.besthits.append(besthit)
                if besthit.score==None:
                    total_score+=childpeak.missing_fragment_score
                    # total_score+=missingfragmentpenalty*weight
                else:
                    total_score+=min(besthit.score,childpeak.missing_fragment_score)
                    # total_score+=min(besthit.score,missingfragmentpenalty)*weight
            hit.score = hit.score + total_score
        return hit

    def add_fragment_data_to_hit(hit):
        if hit.fragment != 0:
            hit.atomstring,hit.atomlist=fragment_engine.get_fragment_info(hit.fragment)
            #except:
            #    exit('failed inchi for: '+atomstring+'--'+str(hit.fragment))
            if len(hit.besthits)>0:
                for childhit in hit.besthits:
                    if childhit != None: # still need to work out how to deal with missed fragments
                        add_fragment_data_to_hit(childhit)

# for peak in self.db_session.query(Peak).filter(Scan.mslevel)==1.filter(Peak.intensity>MSfilter)
    Fragmented=False
    hits=[]
    for peak in peaks:
        if not ((not use_all_peaks) and peak.childscan==None):
            protonation=ionisation_mode-(molformula.find('+')>=0)*1
            deltaH=massmatch(peak,mim,protonation,protonation)
            if type(deltaH)==int:
                if not Fragmented:
                    #sys.stderr.write('\nMetabolite '+str(structure.metid)+': '+str(structure.origin)+' '+str(structure.reactionsequence)+'\n')
                    #sys.stderr.write('Mim: '+str(structure.mim)+'\n')
                    fragment_engine=Fragmentation.FragmentEngine(mol,max_broken_bonds,max_small_losses)
                    #fragment_engine=GrowingEngine(mol)
                    if fragment_engine.accepted():
                        fragment_engine.generate_fragments()
                    #sys.stderr.write('N fragments kept: '+str(len(fragment_engine.fragments))+"\n")
                    Fragmented=True
                if fragment_engine.accepted():
                    hit=gethit(peak,(1<<fragment_engine.get_natoms())-1,0,0,mim,-deltaH)
                    add_fragment_data_to_hit(hit)
                    hits.append(hit)
    return hits



