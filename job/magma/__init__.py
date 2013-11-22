#!/usr/bin/env python

import sys,base64,subprocess,StringIO,time,re,os
import sqlite3,struct,zlib,gzip,copy
import pkg_resources
import numpy
import logging
from lxml import etree
from sqlalchemy import create_engine,and_,desc,distinct
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql import func
from models import Base, Metabolite, Reaction, Scan, Peak, Fragment, Run
import requests,functools,macauthlib #required to update callback url
from requests.auth import AuthBase
import pp
import cPickle as pickle
import types
import pars
from operator import itemgetter

import ConfigParser
config = ConfigParser.ConfigParser()
# default to using rdkit if no config can be found
config.add_section('magma job')
config.set('magma job', 'chemical_engine', 'rdkit')
# read config file from current working directory or users home dir
config.read(['magma_job.ini', os.path.expanduser('~/magma_job.ini')])

if config.get('magma job','chemical_engine')=="rdkit":
    import rdkit_engine as Chem     # Use rdkit_engine
elif config.get('magma job','chemical_engine')=="cdk":
    import cdk_engine               # Use cdk_engine
    Chem=cdk_engine.engine()

"""
RDkit dependencies:
calculate molecular formula
read smiles
generate 2D conformation for those
generate smiles
"""


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
                 max_broken_bonds=3,
                 max_water_losses=1,
                 ms_intensity_cutoff=1e6,
                 msms_intensity_cutoff=5,
                 mz_precision=5.0,
                 mz_precision_abs=0.001,
                 precursor_mz_precision=0.005,
                 use_all_peaks=False,
                 call_back_url=None
                 ):
        return AnnotateEngine(self.db_session,ionisation_mode,skip_fragmentation,max_broken_bonds,max_water_losses,
                 ms_intensity_cutoff,msms_intensity_cutoff,mz_precision,mz_precision_abs,precursor_mz_precision,use_all_peaks,call_back_url)
    def get_select_engine(self):
        return SelectEngine(self.db_session)
    def get_data_analysis_engine(self):
        return DataAnalysisEngine(self.db_session)
    def get_call_back_engine(self,
                 id,
                 key
                 ):
        return CallBackEngine(id,key)
    def commit(self):
        self.db_session.commit()
    def close(self):
        self.db_session.close()

class CallBackEngine(object):
    def __init__(self, url):
        self.access_token=config.get('magma job','macs.id')
        self.mac_key=config.get('magma job','macs.key')
        self.url=url
    def update_callback_url(self,status):
        class HTTPMacAuth(AuthBase):
            """Attaches HTTP Basic Authentication to the given Request object."""
            def __init__(self, id, key):
                self.id = id
                self.key = key
            def __call__(self, r):
                r.headers['Authorization'] = macauthlib.sign_request(r, id=self.id, key=self.key)
                return r

        r = requests.put(self.url, status, auth=HTTPMacAuth(self.access_token, self.mac_key))
        #print r

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

    def add_structure(self,molblock,name,prob,level,sequence,isquery,mim=None,natoms=None,inchikey=None,molform=None,reference=None,logp=None,mass_filter=9999):
        molecule=types.MoleculeType(molblock,name,prob,level,sequence,isquery,mim,natoms,inchikey,molform,reference,logp)
        self.add_molecule(molecule,mass_filter)

    def add_molecule(self,molecule,mass_filter=9999,check_duplicates=True,merge=False):
        if molecule.mim > mass_filter:
            return
        metab=Metabolite(
            mol=unicode(molecule.molblock, 'utf-8', 'xmlcharrefreplace'),
            level=molecule.level,
            probability=molecule.probability,
            reactionsequence=molecule.reactionsequence,
            smiles=molecule.inchikey,
            molformula=molecule.molformula,
            isquery=molecule.isquery,
            origin=unicode(molecule.name, 'utf-8', 'xmlcharrefreplace'),
            nhits=0,
            mim=molecule.mim,
            natoms=molecule.natoms,
            reference=molecule.reference,
            logp=molecule.logp
            )
        if check_duplicates:
            dups=self.db_session.query(Metabolite).filter_by(smiles=molecule.inchikey).all()
            if len(dups)>0:
                if merge:
                    metab=dups[0]
                    metab.origin=unicode(str(metab.origin)+'</br>'+molecule.name, 'utf-8', 'xmlcharrefreplace')
                    metab.reference=molecule.reference
                    metab.probability=molecule.probability
                    metab.reactionsequence=unicode(str(metab.reactionsequence)+'</br>'+molecule.reactionsequence)
                else:
                    if dups[0].probability < molecule.probability:
                        metab.metid=dups[0].metid
                        self.db_session.delete(dups[0])
                        #self.db_session.add(metab)
                        sys.stderr.write('Duplicate structure: - old one removed\n')
                    # TODO remove any fragments related to this structure as well
                    else:
                        sys.stderr.write('Duplicate structure: - kept old one\n')
                        return
        self.db_session.add(metab)
        #sys.stderr.write('Added: '+name+'\n')
        self.db_session.flush()
        return metab.metid

    def metabolize(self,metid,metabolism,endpoints=False):
        # Define if reactions are performed with reactor or with cactvs toolbox
        metabolism_engine=config.get('magma job','metabolism_engine')
        if metabolism_engine=="reactor":
            exec_reactor=pkg_resources.resource_filename( #@UndefinedVariable
                                                      'magma', 'script/reactor')
            exec_reactor+=" -s" # -f 0.15
            if endpoints:
                exec_reactor+=" -m 10"
            else:
                exec_reactor+=" -a -m 1"
            metabolism_files={
                "phase1": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/sygma_rules_GE_0.1.phase1.smirks"),
                "phase2": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/sygma_rules.phase2.smirks"),
                "gut": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/gut.smirks"),
                "digest": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/digest.smirks")
                }
        elif metabolism_engine=="cactvs":
            exec_reactor=pkg_resources.resource_filename( #@UndefinedVariable
                                                      'magma', 'script/csreact')
            if endpoints:
                exec_reactor+=" "+"endpoints"
            else:
                exec_reactor+=" "+"parallel"
            metabolism_files={
                "phase1": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/sygma_rules4.0_GE_0.1.cactvs.phase1.smirks"),
                "phase2": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/sygma_rules4.0.cactvs.phase2_GE0.05.smirks"),
                "gut": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/gut.cactvs.smirks"),
                "digest": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/digest.cactvs.smirks"),
                "peptide": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/peptide.cactvs.smirks"),
                "ptm": pkg_resources.resource_filename( #@UndefinedVariable
                                                           'magma', "data/ptm.cactvs.smirks")
                }
        try:
            parent = self.db_session.query(Metabolite).filter_by(metid=metid).one()
        except:
            print 'Metabolite record ',metid,' does not exist.'
            return
        for m in metabolism.split(','):
            if m in metabolism_files:
                if metabolism_engine=="reactor":
                    exec_reactor+=" -q"
                exec_reactor+=" "+metabolism_files[m]
        sys.stderr.write(exec_reactor + '\n')

        reactor=subprocess.Popen(exec_reactor, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        reactor.stdin.write(parent.mol+'$$$$\n')
        reactor.stdin.close()

        metids=set()
        line=reactor.stdout.readline()
        while line != "":
            splitline=line[:-1].split(" {>")
            if len(splitline) == 1:
                reaction='PARENT'
                name=line[:-1]
            else:
                reaction=" + ".join([r[:-1] for r in splitline[1:]])
                name=splitline[0]
            mol=name+'\n'
            isquery=0
            while line != 'M  END\n':
                line=reactor.stdout.readline()
                mol+=line
            line=reactor.stdout.readline()
            while line != '$$$$\n' and line != "":
                #if line=='> <Probability>\n':
                #    prob=float(reactor.stdout.readline())
                #elif line=='> <Level>\n':
                #    level=int(reactor.stdout.readline())
                if line=='> <ReactionSequence>\n':
                    line=reactor.stdout.readline()
                    reaction=line[:-1]
                    line=reactor.stdout.readline()
                    while line != '\n':
                        reaction+=' + '+line[:-1]
                        line=reactor.stdout.readline()
                line=reactor.stdout.readline()
            if reaction!='PARENT':
                molecule=types.MoleculeType(mol,name,0,0,str(metid)+'&gt&gt'+reaction,isquery)
                new_metid=self.add_molecule(molecule,merge=True)
                metids.add(new_metid)
                reactid=self.db_session.query(Reaction.reactid).filter(Reaction.reactant==metid,Reaction.product==new_metid,Reaction.name==reaction).all()
                if len(reactid)==0:
                    parent.reactionsequence+=unicode('</br>'+str(new_metid)+'&lt&lt'+reaction)
                    react=Reaction(
                        reactant=metid,
                        product=new_metid,
                        name=reaction
                        )
                    self.db_session.add(react)
                    self.db_session.flush()
            elif endpoints:
                metids.add(metid)
            line=reactor.stdout.readline()
        self.db_session.add(parent)
        reactor.stdout.close()
        self.db_session.commit()
        if len(metids)==0: # this might be the case with cactvs engine
            metids.add(metid)
        return metids

    def metabolize_all(self,metabolism,endpoints=False):
        logging.warn('Metabolize all')
        parentids = self.db_session.query(Metabolite.metid).all()
        # print parentids
        metids=set([])
        for parentid, in parentids:
            metids|=self.metabolize(parentid,metabolism,endpoints)
        return metids

    def run_scenario(self, scenario):
        print scenario
        metids=set()
        old_metids=set()
        for step in range(len(scenario)):
            input,action,value = scenario[step]
            print "----- Scenario, step ",step,"-----"
            endpoints=False
            if input=="previous":
                metids |= old_metids
            old_metids=metids
            if action=='mass_filter':
                if input=='all':
                    result=self.db_session.query(Metabolite.metid).filter(Metabolite.mim<float(value)).all()
                else:
                    result=self.db_session.query(Metabolite.metid).filter(Metabolite.mim<float(value),Metabolite.metid.in_(metids)).all()
                    print "from ",len(metids),"processed compounds"
                metids={x[0] for x in result} #set comprehension
                print len(metids),"compounds were selected with mass <",value
            else:
                if value=='end':
                    endpoints=True
                    value=1
                if input=='all':
                    metids=set()
                    new_metids = self.metabolize_all(action,endpoints)
                    print "All metabolites",
                else:
                    print len(metids),"metabolites",
                    new_metids=set()
                    for metid in metids:
                        new_metids |= self.metabolize(metid,action,endpoints)
                active_metids=new_metids.difference(metids)
                metids=new_metids
                for it in range(1,int(value)):
                    new_metids=set()
                    for metid in active_metids:
                        new_metids |= self.metabolize(metid,action,endpoints)
                    active_metids=new_metids.difference(metids)
                    metids |= new_metids
                print 'were metabolized according to',action,'rules'
                print '"Active" metabolites:',len(metids)
            # print metids
            print ""
 
    def retrieve_structures(self,mass):
        dbfilename = '/home/ridderl/chebi/ChEBI_complete_3star.sqlite'
        conn = sqlite3.connect(dbfilename)
        c = conn.cursor()
        result = c.execute('SELECT * FROM molecules WHERE mim BETWEEN ? AND ?' , (mass-0.01,mass+0.01))
        for (id,mim,molblock,smiles,chebi_name) in result:
            self.add_molecule(zlib.decompress(molblock),str(chebi_name),1.0,1,"",1)

    def get_metids_with_mass_less_then(self,mass,metids=None):
        if metids==None:
            result=self.db_session.query(Metabolite.metid).filter(Metabolite.mim<mass).all()
        else:
            result=self.db_session.query(Metabolite.metid).filter(Metabolite.mim<mass,Metabolite.metid.in_(metids)).all()
        metids={x[0] for x in result} #set comprehension
        return metids

class MsDataEngine(object):
    def __init__(self,db_session,abs_peak_cutoff,rel_peak_cutoff,max_ms_level):
        self.db_session = db_session
        try:
            rundata=self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.abs_peak_cutoff == None:
            rundata.abs_peak_cutoff=abs_peak_cutoff
#        if rundata.rel_peak_cutoff == None:
#            rundata.rel_peak_cutoff=rel_peak_cutoff
        if rundata.max_ms_level == None:
            rundata.max_ms_level=max_ms_level
        self.db_session.add(rundata)
        self.db_session.commit()
        self.abs_peak_cutoff=rundata.abs_peak_cutoff
#        self.rel_peak_cutoff=rundata.rel_peak_cutoff
        self.max_ms_level=rundata.max_ms_level

    def store_mzxml_file(self,mzxml_file,scan_filter=None):
        logging.warn('Store mzxml file')
        rundata=self.db_session.query(Run).one()
        if rundata.ms_filename != None:
            sys.exit('Attempt to read MS data twice')
        rundata.ms_filename=mzxml_file
        self.db_session.add(rundata)
        self.ms_filename=mzxml_file
        tree=etree.parse(mzxml_file)
        root=tree.getroot()
        namespace='{'+root.nsmap[None]+'}'
        mzxml_query=namespace+"msRun/"+namespace+"scan"
        prec_scans=[] # in case of non-hierarchical mzXML, also find child scans
        for mzxmlScan in root.findall(mzxml_query):
            try:
                if scan_filter == None or mzxmlScan.attrib['num'] == scan_filter or mzxmlScan.find(namespace+'precursorMz').attrib['precursorScanNum'] in prec_scans:
                    self.store_mzxml_scan(mzxmlScan,0,namespace)
                    prec_scans.append(mzxmlScan.attrib['num']) # in case of non-hierarchical mzXML, also find child scans
            except:
                pass
        self.db_session.commit()

    def store_mzxml_scan(self,mzxmlScan,precScan,namespace):
        if mzxmlScan.attrib['peaksCount']=='0':
            return
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
                decoded = base64.decodestring(child.text)
                try:
                    if child.attrib['compressionType']=='zlib':
                        decoded=zlib.decompress(decoded)
                except:
                    pass
                self.store_mzxml_peaks(scan,decoded)
            if child.tag == namespace+'scan' and int(child.attrib['msLevel'])<=self.max_ms_level:
#                if scan.mslevel == 1:
#                    cutoff=0.0
#                else:
#                    cutoff=max(scan.basepeakintensity*self.rel_peak_cutoff,self.abs_peak_cutoff)
#                if float(child.attrib['basePeakIntensity'])>cutoff:
                self.store_mzxml_scan(child,scan.scanid,namespace)
        self.db_session.add(scan)
        self.db_session.flush()

    def store_mzxml_peaks(self,scan,decoded):
#        dbpeak_cutoff=scan.basepeakintensity*self.rel_peak_cutoff
#        if dbpeak_cutoff<self.abs_peak_cutoff:
#            dbpeak_cutoff=self.abs_peak_cutoff
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

    def store_manual_tree(self,manual_tree,tree_type):
        # tree_type: 0 for mass tree, -1 and 1 for formula trees with negative or positive ionisation mode respectively
        tree_string=''.join(open(manual_tree).read().split()) # remove whitespaces (' ','\t','\n',etc) from tree_string
        tree_list=re.split('([\,\(\)])',tree_string)
        self.global_scanid = 1
        self.store_manual_subtree(tree_list,0,0,0,1,tree_type)

    def store_manual_subtree(self,tree_list,precursor_scanid,precursor_mz,precursor_intensity,mslevel,tree_type):
        lowmz=None
        highmz=None
        basepeakmz=None
        basepeakintensity=None
        scanid=self.global_scanid
        npeaks=0
        while len(tree_list)>0 and tree_list[0]!=')':
            tree_item=tree_list.pop(0)
            if tree_item.find(':')>=0:
                mz,intensity=tree_item.split(':')
                if tree_type != 0:
                    mz=self.mass_from_formula(mz)-tree_type*pars.elmass
                self.db_session.add(Peak(scanid=scanid,mz=mz,intensity=intensity))
                npeaks+=1
                if lowmz==None or mz<lowmz:
                    lowmz=mz
                if highmz==None or mz>highmz:
                    highmz=mz
                if basepeakintensity==None or intensity>basepeakintensity:
                    basepeakmz=mz
                    basepeakintensity=intensity
            elif tree_item=='(':
                self.global_scanid+=1
                self.store_manual_subtree(tree_list,scanid,mz,intensity,mslevel+1,tree_type)
            elif tree_item!=',' and tree_item!='':
                exit('Corrupt Tree format ...')
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

    def mass_from_formula(self,form):
        mass=0.0
        while len(form)>0:
            if form[:2] in pars.mims:
                m=pars.mims[form[:2]]
                form=form[2:]
            elif form[:1] in pars.mims:
                m=pars.mims[form[:1]]
                form=form[1:]
            else:
                exit('Element not allowed in formula tree: '+form)
            x=0
            while len(form)>x and form[x] in '0123456789':
                x+=1
            if x>0:
                n=int(form[:x])
            else:
                n=1
            mass+=m*n
            form=form[x:]
        return mass

class AnnotateEngine(object):
    def __init__(self,db_session,ionisation_mode,skip_fragmentation,max_broken_bonds,max_water_losses,
                 ms_intensity_cutoff,msms_intensity_cutoff,mz_precision,mz_precision_abs,
                 precursor_mz_precision,use_all_peaks,call_back_url=None):
        self.db_session = db_session
        mz_precision_abs=max(mz_precision_abs,0.000001)
        precursor_mz_precision=max(precursor_mz_precision,0.000001)
        # a small mz_precision_abs is required, even when matching theoretical masses, because of finite floating point precision
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
        if rundata.max_water_losses == None:
            rundata.max_water_losses=max_water_losses
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
        self.max_water_losses=rundata.max_water_losses
        self.ms_intensity_cutoff=rundata.ms_intensity_cutoff
        self.msms_intensity_cutoff=rundata.msms_intensity_cutoff
        self.mz_precision=rundata.mz_precision
        self.precision=1+rundata.mz_precision/1e6
        self.mz_precision_abs=rundata.mz_precision_abs
        self.precursor_mz_precision=rundata.precursor_mz_precision
        self.use_all_peaks=rundata.use_all_peaks

        self.scans=[]

        if call_back_url != None:
            self.call_back_engine=CallBackEngine(call_back_url)
        else:
            self.call_back_engine=None

    def build_spectrum(self,dbscan):
        scan=types.ScanType(dbscan.scanid,dbscan.mslevel)
        if scan.mslevel==1:
            cutoff=self.ms_intensity_cutoff
        else:
            cutoff=dbscan.basepeakintensity*self.msms_intensity_cutoff/100
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

    def get_db_candidates(self,query_engine,max_mim=""):
        if max_mim=='':
            max_mim='1200'
        mmim=int(float(max_mim)*1e6)

        struct_engine = StructureEngine(self.db_session,"",0)

        # build sorted list of query masses
        mzs=[]
        for scan in self.scans:
            # if self.db_session.query(Fragment.fragid).filter(Fragment.scanid==scan.scanid).count() > 0: # tijdelijk!
            for peak in scan.peaks:
                if not ((not self.use_all_peaks) and peak.childscan==None):
                    mzs.append(peak.mz)
        mzs.sort()
        # build non-overlapping set of queries around these masses
        queries=[[0,0]]
        cqueries=[[0,0]] # queries of charged molecules
        for mz in mzs:
            ql=int(1e6*(min(mz/self.precision,mz-self.mz_precision_abs)-self.ionisation_mode*(pars.Hmass-pars.elmass)))
            qh=int(1e6*(max(mz*self.precision,mz+self.mz_precision_abs)-self.ionisation_mode*(pars.Hmass-pars.elmass)))
            #ql=int(1e6*(mz/self.precision-self.ionisation_mode*(pars.Hmass-pars.elmass))) # tijdelijk for thermo
            #qh=int(1e6*(mz*self.precision-self.ionisation_mode*(pars.Hmass-pars.elmass))) # tijdelijk for thermo
            if ql < mmim:
                if qh > mmim:
                    qh == mmim
                if queries[-1][0] <= ql <= queries[-1][1]:
                    queries[-1][1]=qh
                else:
                    queries.append([ql,qh])
            ql=int(1e6*(min(mz/self.precision,mz-self.mz_precision_abs)+self.ionisation_mode*pars.elmass))
            qh=int(1e6*(max(mz*self.precision,mz+self.mz_precision_abs)+self.ionisation_mode*pars.elmass))
            if ql < mmim:
                if qh > mmim:
                    qh == mmim
                if cqueries[-1][0] <= ql <= cqueries[-1][1]:
                    cqueries[-1][1]=qh
                else:
                    cqueries.append([ql,qh])

        # in case of an empty database, no check for existing duplicates needed
        check_duplicates = (self.db_session.query(Metabolite.metid).count() > 0)
        print 'check: ',check_duplicates

        # All candidates are stored in dbsession, resulting metids are returned
        metids=set([])
        for ql,qh in queries:
            result=query_engine.query_on_mim(ql,qh,0)
            for molecule in result:
                metid=struct_engine.add_molecule(molecule,check_duplicates=check_duplicates,merge=True)
                metids.add(metid)
            print str(ql)+','+str(qh)+' --> '+str(len(metids))+' candidates'
        for ql,qh in cqueries:
            result=query_engine.query_on_mim(ql,qh,self.ionisation_mode)
            for molecule in result:
                metid=struct_engine.add_molecule(molecule,check_duplicates=check_duplicates,merge=True)
                metids.add(metid)
            print str(ql)+','+str(qh)+' --> '+str(len(metids))+' candidates'
        self.db_session.commit()
        return metids

    def search_structures(self,metids=None,ncpus=1,fast=False,time_limit=None):
        logging.warn('Searching structures')
        global fragid
        fragid=self.db_session.query(func.max(Fragment.fragid)).scalar()
        if fragid == None:
            fragid = 0
        ppservers = ()
        logging.warn('calculating on '+str(ncpus)+' cpus !!!')
        job_server = pp.Server(ncpus, ppservers=ppservers)
        if metids==None:
            if time_limit == None:
                metabdata=self.db_session.query(Metabolite.metid).order_by(desc(Metabolite.metid)).all()
            else:
                metabdata=self.db_session.query(Metabolite.metid).order_by(Metabolite.probability).all()
            metids=[x[0] for x in metabdata]
        # annotate metids in chunks of 500 to avoid errors in db_session.query and memory problems during parallel processing
        total_frags=0
        total_metids = len(metids)
        start_time=time.time()
        update_time=start_time
        update_interval=1 #send update to call_back_url every second
        while len(metids)>0:
            ids=set([])
            while len(ids)<500 and len(metids)>0:
                ids.add(metids.pop())
            structures = self.db_session.query(Metabolite).filter(Metabolite.metid.in_(ids)).all()
            jobs=[]
            for structure in structures:
                # skip molecule if it has already been used for annotation
                if self.db_session.query(Fragment.fragid).filter(Fragment.metid == structure.metid).count() > 0:
                    continue
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
                    if fast and structure.natoms<=64:
                        fragmentation_module='magma.fragmentation_cy'
                    else:
                        fragmentation_module='magma.fragmentation_py'
                    jobs.append((structure,
                       job_server.submit(search_structure,(structure.mol,
                                  structure.mim,
                                  structure.molformula,
                                  peaks,
                                  self.max_broken_bonds,
                                  self.max_water_losses,
                                  self.precision,
                                  self.mz_precision_abs,
                                  self.use_all_peaks,
                                  self.ionisation_mode,
                                  self.skip_fragmentation,
                                  (fast and structure.natoms<=64),
                                  config.get('magma job','chemical_engine')
                                  ),(),(
                                  "magma.types",
                                  "magma.pars",
                                  fragmentation_module
                                  )
                               )))
            count=0
            for structure,job in jobs:
                sys.stderr.write('Metabolite '+str(structure.metid))
                raw_result=job(raw_result=True)
                result,sout = pickle.loads(raw_result)
                print sout
                (hits,frags)=result
                total_frags+=frags
                sys.stderr.write(' -> '+str(frags)+' fragments: '+str(structure.origin.encode('utf8'))+'\n')
                structure.nhits=len(hits)
                self.db_session.add(structure)
                for hit in hits:
                    sys.stderr.write('Scan: '+str(hit.scan)+' - Mz: '+str(hit.mz)+' - ')
                    sys.stderr.write('Score: '+str(hit.score)+'\n')
                    self.store_hit(hit,structure.metid,0)
                self.db_session.flush()
                count+=1
                elapsed_time=time.time()-start_time
                if self.call_back_engine != None:
                    update_last=time.time()-update_time
                    if update_last > update_interval: # update status every second
                        update_string='<h1>Progress: %d / %d candidate molecules processed  (%d%%)</h1>' % (total_metids-len(metids)-len(ids)+count,total_metids,100.0*(total_metids-len(metids)-len(ids)+count)/total_metids)
                        update_string+='<h3>Time: %02d:%02d:%02d' % (elapsed_time//3600,(elapsed_time%3600)//60,elapsed_time%60)
                        if time_limit != None:
                            update_string+=' / max. %02d:%02d:00 (%d%%)' % (time_limit//60,time_limit%60,elapsed_time/time_limit/60*100)
                        update_string+='</h3>'
                        self.call_back_engine.update_callback_url(update_string)
                        update_time = update_time + update_last//update_interval*update_interval
                if time_limit and elapsed_time > time_limit * 60:
                    metids=[] # break out of while-loop
                    sys.stderr.write('\nThe time limit of '+str(time_limit)+' minutes has been exceeded!\n\n')
                    break
            self.db_session.commit()
        print total_frags,'fragments in total.'

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
            deltappm=deltappm,
            formula=hit.formula
            ))
        if len(hit.besthits)>0:
            for childhit in hit.besthits:
                if childhit != None: # still need to work out how to deal with missed fragments
                    self.store_hit(childhit,metid,currentFragid)

class PubChemEngine(object):
    def __init__(self,dbfilename='',max_64atoms=False,min_refscore=''):
        if dbfilename=='':
            dbfilename=config.get('magma job','structure_database.pubchem')
        self.where=''
        if min_refscore!='':
            self.where += ' AND refscore >= '+min_refscore
        if max_64atoms==True:
            self.where += ' AND natoms <= 64'
        self.conn = sqlite3.connect(dbfilename)
        self.conn.text_factory=str
        self.c = self.conn.cursor()
    def query_on_mim(self,low,high,charge):
        result=self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge,low,high))
        molecules=[]
        for (cid,mim,charge,natoms,molblock,inchikey,molform,name,refscore,logp) in result:
            molecules.append(types.MoleculeType(
                           molblock=zlib.decompress(molblock),
                           name=name+' ('+str(cid)+')',
                           mim=float(mim/1e6),
                           natoms=natoms,
                           molform=molform,
                           inchikey=inchikey,
                           prob=refscore,
                           level=1,
                           sequence="",
                           isquery=1,
                           reference='<a href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&cmd=Link&LinkName=pccompound_pccompound_sameisotopic_pulldown&from_uid='+\
                                     str(cid)+'" target="_blank">'+str(cid)+' (PubChem)</a>',
                           logp=float(logp)/10.0,
                           ))
        return molecules

class KeggEngine(object):
    def __init__(self,dbfilename='',max_64atoms=False):
        if dbfilename=='':
            dbfilename=config.get('magma job','structure_database.kegg')
        self.where=''
        if max_64atoms==True:
            self.where += ' AND natoms <= 64'
        self.conn = sqlite3.connect(dbfilename)
        self.conn.text_factory=str
        self.c = self.conn.cursor()
    def query_on_mim(self,low,high,charge):
        result=self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge,low,high))
        molecules=[]
        for (cid,mim,charge,natoms,molblock,inchikey,molform,name,reference,logp) in result:
            keggids=reference.split(',')
            keggrefs='<a href="http://www.genome.jp/dbget-bin/www_bget?cpd:'+keggids[0]+'" target="_blank">'+keggids[0]+' (Kegg)</a>'
            for keggid in keggids[1:]:
                keggrefs+='<br><a href="http://www.genome.jp/dbget-bin/www_bget?cpd:'+keggid+'" target="_blank">'+keggid+' (Kegg)</a>'
            molecules.append(types.MoleculeType(
                           molblock=zlib.decompress(molblock),
                           name=name+' ('+str(cid)+')',
                           mim=float(mim/1e6),
                           natoms=natoms,
                           molform=molform,
                           inchikey=inchikey,
                           prob=None,
                           level=1,
                           sequence="",
                           isquery=1,
                           reference=keggrefs,
                           logp=float(logp)/10.0,
                           ))
        return molecules

class HmdbEngine(object):
    def __init__(self,dbfilename='',max_64atoms=False):
        if dbfilename=='':
            dbfilename=config.get('magma job','structure_database.hmdb')
        self.where=''
        if max_64atoms==True:
            self.where += ' AND natoms <= 64'
        self.conn = sqlite3.connect(dbfilename)
        self.conn.text_factory=str
        self.c = self.conn.cursor()
    def query_on_mim(self,low,high,charge):
        result=self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge,low,high))
        molecules=[]
        for (cid,mim,charge,natoms,molblock,inchikey,molform,name,reference,logp) in result:
            hmdb_ids=reference.split(',')
            hmdb_refs='<a href="http://www.hmdb.ca/metabolites/'+hmdb_ids[0]+'" target="_blank">'+hmdb_ids[0]+' (HMDB)</a>'
            for hmdb_id in hmdb_ids[1:]:
                hmdb_refs+='<br><a href="http://www.hmdb.ca/metabolites/'+hmdb_id+'" target="_blank">'+hmdb_id+' (HMDB)</a>'
            molecules.append(types.MoleculeType(
                           molblock=zlib.decompress(molblock),
                           name=name+' ('+str(cid)+')',
                           mim=float(mim/1e6),
                           natoms=natoms,
                           molform=molform,
                           inchikey=inchikey,
                           prob=None,
                           level=1,
                           sequence="",
                           isquery=1,
                           reference=hmdb_refs,
                           logp=float(logp)/10.0,
                           ))
        return molecules

class MetaCycEngine(object):
    def __init__(self,dbfilename='',max_64atoms=False):
        if dbfilename=='':
            dbfilename=config.get('magma job','structure_database.metacyc')
        self.where=''
        if max_64atoms==True:
            self.where += ' AND natoms <= 64'
        self.conn = sqlite3.connect(dbfilename)
        self.conn.text_factory=str
        self.c = self.conn.cursor()
    def query_on_mim(self,low,high,charge):
        result=self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge,low,high))
        molecules=[]
        for (cid,mim,charge,natoms,molblock,inchikey,molform,name,reference,logp) in result:
            metacyc_ids=reference.split(',')
            metacyc_refs='<a href="http://www.biocyc.org/META/NEW-IMAGE?type=COMPOUND&object='+metacyc_ids[0]+'" target="_blank">'+metacyc_ids[0]+' (MetaCyc)</a>'
            for metacyc_id in metacyc_ids[1:]:
                metacyc_refs+='<br><a href="http://www.biocyc.org/META/NEW-IMAGE?type=COMPOUND&object='+metacyc_id+'" target="_blank">'+metacyc_id+' (MetaCyc)</a>'
            molecules.append(types.MoleculeType(
                           molblock=zlib.decompress(molblock),
                           name=name,
                           mim=float(mim/1e6),
                           natoms=natoms,
                           molform=molform,
                           inchikey=inchikey,
                           prob=None,
                           level=1,
                           sequence="",
                           isquery=1,
                           reference=metacyc_refs,
                           logp=float(logp)/10.0,
                           ))
        return molecules

class SelectEngine(object):
    def __init__(self,db_session):
        self.db_session = db_session
        mz_precision,mz_precision_abs=self.db_session.query(Run.mz_precision,Run.mz_precision_abs).all()[0]
        self.precision=1+mz_precision/1e6
        self.mz_precision_abs=mz_precision_abs

    def select_fragment(self,fragid):
        child_frags = self.db_session.query(Fragment.fragid).filter(Fragment.parentfragid==fragid).all()
        if len(child_frags) == 0:
            exit('Fragment not fragmented')
        fragmented_fragids = self.db_session.query(distinct(Fragment.parentfragid))
        smiles = self.db_session.query(Fragment.inchikey).filter(Fragment.fragid==fragid).all()
        sys.stderr.write(str(smiles)+'\n')
        metids = self.db_session.query(Fragment.metid).filter(Fragment.inchikey==smiles[0][0], Fragment.fragid.in_(fragmented_fragids))
        # print frags.all()
        self.db_session.query(Fragment).filter(~Fragment.metid.in_(metids)).delete(synchronize_session='fetch')
        self.db_session.query(Metabolite).filter(~Metabolite.metid.in_(metids)).delete(synchronize_session='fetch')
        self.db_session.commit()
        ref_scan = self.db_session.query(distinct(Fragment.scanid)).filter(Fragment.parentfragid==fragid).all()[0][0]
        fragids = self.db_session.query(Fragment.fragid).filter(Fragment.inchikey==smiles[0][0], Fragment.fragid.in_(fragmented_fragids))
        query_scans = self.db_session.query(distinct(Fragment.scanid)).filter(Fragment.parentfragid.in_(fragids)).all()
        for query_scan, in query_scans:
            dot_product=self.dot_product_scans(ref_scan,query_scan)
            compounds = self.db_session.query(Metabolite).filter(Metabolite.metid.in_(
                            self.db_session.query(Fragment.metid).filter(Fragment.scanid==query_scan))).all()
            for compound in compounds:
                compound.reactionsequence+='Scan: '+str(query_scan)+' - Similarity: '+str(dot_product)+'\n'
            self.db_session.add(compound)
        self.db_session.commit()

    def dot_product_scans(self,ref_scan,query_scan):
        ref_peaks=self.db_session.query(Peak.mz,Peak.intensity).filter(Peak.scanid==ref_scan).all()
        query_peaks=self.db_session.query(Peak.mz,Peak.intensity).filter(Peak.scanid==query_scan).all()
        ref_peaks_sorted=sorted(ref_peaks, key=itemgetter(1), reverse=True)        # Start peak matching with most intense peaks first
        query_peaks_sorted=sorted(query_peaks, key=itemgetter(1), reverse=True)
        dot_product=None
        if len(query_peaks) > 0:
            matched_intensities=[]
            for pr in ref_peaks_sorted:
                #print ref_peaks_sorted
                #print '-------------'
                #print query_peaks_sorted
                #print '============='
                ll = min(pr[0]/self.precision,pr[0]-self.mz_precision_abs)
                hl = max(pr[0]*self.precision,pr[0]+self.mz_precision_abs)
                min_delta_intensity=None
                match=None
                for pq in query_peaks_sorted:
                    if ll < pq[0] < hl:
                        delta_intensity=abs(pq[1]-pr[1])
                        if min_delta_intensity==None or min_delta_intensity>delta_intensity:
                            min_delta_intensity=delta_intensity
                            match=pq
                if match == None:
                    matched_intensities.append((pr[1],0))
                else:
                    matched_intensities.append((pr[1],match[1]))
                    query_peaks_sorted.remove(match)
            for pq in query_peaks_sorted:
                matched_intensities.append((0,pq[1]))
            #print matched_intensities
            srq=0
            sr=0
            sq=0
            for match in matched_intensities:
                srq+=(match[0]*match[1])**0.5
                sr+=match[0]
                sq+=match[1]
            dot_product=srq**2/(sr*sq)
        return dot_product
        

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

    def write_SDF(self,file=sys.stdout,molecules=None,columns=None,sortcolumn=None,descend=False):
        if molecules==None:
            if descend:
                molecules=self.db_session.query(Metabolite).order_by(desc(sortcolumn)).all()
            else:
                molecules=self.db_session.query(Metabolite).order_by(sortcolumn).all()
        for molecule in molecules:
            file.write(molecule.mol)
            if columns==None:
                columns=dir(molecule)
            for column in columns:
                if column[:1] != '_' and column != 'mol' and column != 'metadata':
                    file.write('> <'+column+'>\n'+str(molecule.__getattribute__(column))+'\n\n')
            file.write('$$$$\n')

    def write_network1(self,filename):
        f=open(filename+'.sif','w')
        assigned_metids=self.db_session.query(distinct(Peak.assigned_metid))
        written_metids=set([])
        for reactant,product in self.db_session.query(Reaction.reactant,Reaction.product).filter(Reaction.reactant.in_(assigned_metids) | Reaction.product.in_(assigned_metids)).all():
            f.write(str(reactant)+' pp '+str(product)+'\n')
            written_metids.add(reactant)
            written_metids.add(product)
        f.close()
        f=open(filename+'.txt','w')
        assigned=assigned_metids.all()
        for metid in written_metids:
            if (metid,) in assigned:
                f.write(str(metid)+' assigned\n')
            else:
                f.write(str(metid)+' unassigned\n')

    def write_network2(self,filename):
        nodes={}
        for reactant,product in self.db_session.query(Reaction.reactant,Reaction.product).all():
            r=int(reactant)
            p=int(product)
            if r not in nodes:
                nodes[r]=set([])
            if p not in nodes:
                nodes[p]=set([])
            nodes[r].add(p)
            nodes[p].add(r)
        result=self.db_session.query(distinct(Peak.assigned_metid)).all()
        assigned_metids=[x[0] for x in result]
        start_compound={}
        result=self.db_session.query(Metabolite.metid,Metabolite.isquery).all()
        for metid,isquery in result:
            start_compound[metid]=isquery
        print result
        print assigned_metids
        nnodes=len(nodes)+1
        while len(nodes) < nnodes:
            nnodes = len(nodes)
            print nnodes
            nodekeys=nodes.keys()
            for n in nodekeys:
                if n not in assigned_metids and start_compound[n] == False:
                    if len(nodes[n]) == 1: # and list(nodes[n])[0] not in assigned_metids:
                        nodes[list(nodes[n])[0]].remove(n)
                        del nodes[n]
                    else:
                        if len(nodes[n]) == 2:
                            tmpnode=list(nodes[n])
                            #if len(nodes[tmpnode[0]] & nodes[tmpnode[1]]) > 1 and tmpnode[0] not in assigned_metids and tmpnode[1] not in assigned_metids:
                            if len(nodes[tmpnode[0]] & nodes[tmpnode[1]]) > 1 or \
                                    len(nodes[tmpnode[0]] & nodes[n]) > 0 or \
                                    len(nodes[tmpnode[1]] & nodes[n]) > 0:
                                for c in nodes[n]:
                                    nodes[c].remove(n)
                                del nodes[n]
        f=open(filename+'.sif','w')
        written_metids=set([])
        connections=[]
        for n in nodes:
            for c in nodes[n]:
                l = set([n,c])
                if l not in connections:
                    f.write(str(n)+' pp '+str(c)+'\n')
                    connections.append(l)
            written_metids.add(n)
        f.close()
        f=open(filename+'.txt','w')
        for metid in written_metids:
            if (metid) in assigned_metids:
                f.write(str(metid)+" "+str(start_compound[metid]+2)+'\n')
            else:
                f.write(str(metid)+" "+str(start_compound[metid]+0)+'\n')

def search_structure(mol,mim,molformula,peaks,max_broken_bonds,max_water_losses,precision,mz_precision_abs,use_all_peaks,ionisation_mode,skip_fragmentation,fast,chem_engine):
    pars=magma.pars
    if fast:
        Fragmentation=magma.fragmentation_cy
    else:
        Fragmentation=magma.fragmentation_py

    def massmatch(peak,mim,protonation):
        lowmz=min(peak.mz/precision,peak.mz-mz_precision_abs)
        highmz=max(peak.mz*precision,peak.mz+mz_precision_abs)
        #lowmz=peak.mz/precision #for peptides file
        #highmz=peak.mz*precision #for peptides file
        return lowmz <= (mim+protonation*pars.Hmass-ionisation_mode*pars.elmass) <= highmz

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
                if besthit.score==None:
                    total_score+=childpeak.missing_fragment_score
                    # total_score+=missingfragmentpenalty*weight
                else:
                    hit.besthits.append(besthit)
                    total_score+=min(besthit.score,childpeak.missing_fragment_score)
                    # total_score+=min(besthit.score,missingfragmentpenalty)*weight
            hit.score = hit.score + total_score
        return hit

    def add_fragment_data_to_hit(hit):
        if hit.fragment != 0:
            hit.atomstring,hit.atomlist,hit.formula,hit.inchikey=fragment_engine.get_fragment_info(hit.fragment,hit.deltaH)
            #except:
            #    exit('failed inchi for: '+atomstring+'--'+str(hit.fragment))
            if len(hit.besthits)>0:
                for childhit in hit.besthits:
                    if childhit != None: # still need to work out how to deal with missed fragments
                        add_fragment_data_to_hit(childhit)

# for peak in self.db_session.query(Peak).filter(Scan.mslevel)==1.filter(Peak.intensity>MSfilter)
    Fragmented=False
    hits=[]
    frags=0
    charge=0
    charge+=-1*(molformula[-1]=='-')+1*(molformula[-1]=='+') # derive charge from molecular formula
    for peak in peaks:
        if not ((not use_all_peaks) and peak.childscan==None):
            protonation=ionisation_mode-charge
            if massmatch(peak,mim,protonation):
                if not Fragmented:
                    #sys.stderr.write('\nMetabolite '+str(structure.metid)+': '+str(structure.origin)+' '+str(structure.reactionsequence)+'\n')
                    #sys.stderr.write('Mim: '+str(structure.mim)+'\n')
                    fragment_engine=Fragmentation.FragmentEngine(mol,max_broken_bonds,max_water_losses,ionisation_mode,skip_fragmentation)
                    #fragment_engine=GrowingEngine(mol)
                    if fragment_engine.accepted():
                        frags=fragment_engine.generate_fragments()
                    #sys.stderr.write('N fragments kept: '+str(len(fragment_engine.fragments))+"\n")
                    Fragmented=True
                if fragment_engine.accepted():
                    hit=gethit(peak,(1<<fragment_engine.get_natoms())-1,0,0,mim,-protonation)
                    add_fragment_data_to_hit(hit)
                    hits.append(hit)
    return (hits,frags)



