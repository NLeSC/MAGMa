#!/usr/bin/env python

import sys,base64,subprocess,StringIO,time
import sqlite3,struct,zlib,gzip,copy
import pkg_resources
import numpy
import logging
# from rdkit import Chem, Geometry
# from rdkit.Chem import AllChem, Descriptors
from lxml import etree
from sqlalchemy import create_engine,and_,desc
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql import func
from models import Base, Metabolite, Scan, Peak, Fragment, Run
import jpype,glob,os
import pp
import cPickle as pickle
import types

"""
RDkit dependencies:
calculate molecular formula
read smiles
generate 2D conformation for those
generate smiles
"""

#jars= glob.glob('/home/ridderl/cdk/Marijn_jars/*jar')
jars= ('/home/ridderl/cdk/cdk-1.4.13.jar',)
classpath = ":".join([ os.path.abspath(jar) for jar in jars])
os.environ['JAVA_HOME'] = '/usr/lib/jvm/java-6-openjdk'
jpype.startJVM(jpype.getDefaultJVMPath(),"-ea", "-Djava.class.path="+classpath)
#
#cdk = jpype.JPackage("org").openscience.cdk
#java = jpype.java


class CDKengine(object):
    def __init__(self):
#        jars= ('/home/ridderl/cdk/cdk-1.4.13.jar',)
#        classpath = ":".join([ os.path.abspath(jar) for jar in jars])
#        os.environ['JAVA_HOME'] = '/usr/lib/jvm/java-6-openjdk'
#        jpype.startJVM(jpype.getDefaultJVMPath(),"-ea", "-Djava.class.path="+classpath)
        self.cdk = jpype.JPackage("org").openscience.cdk
        self.java = jpype.java

        self.builder = self.cdk.DefaultChemObjectBuilder.getInstance()
        self.sp = self.cdk.smiles.SmilesParser(self.builder)
        self.sp.setPreservingAromaticity(True)
        self.sg = self.cdk.smiles.SmilesGenerator()
        self.sg.setUseAromaticityFlag(True)
        self.isof=self.cdk.config.IsotopeFactory.getInstance(self.builder)
        self.Hmass=self.isof.getMajorIsotope('H').getExactMass().floatValue()
        self.acm = self.cdk.tools.manipulator.AtomContainerManipulator
    def MolToMolBlock(self,molecule):
        mol2stringio = self.java.io.StringWriter()
        mol2writer = self.cdk.io.MDLV2000Writer(mol2stringio)
        try:
            dbst = self.cdk.smiles.DeduceBondSystemTool()
            molecule = dbst.fixAromaticBondOrders(molecule)
        except:
            pass
        mol2writer.write(molecule)
        # mol2writer.close()
        return mol2stringio.toString()
    def MolFromMolBlock(self,mol_block):
        stringio2mol = self.java.io.StringReader(mol_block)
        reader2mol = self.cdk.io.MDLReader(stringio2mol)
        molecule=self.cdk.Molecule()
        reader2mol.read(molecule)
        self.cdk.tools.manipulator.AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule)
        #self.cdk.aromaticity.CDKHueckelAromaticityDetector.detectAromaticity(molecule) #
        ha=self.cdk.tools.CDKHydrogenAdder.getInstance(self.builder)
        ha.addImplicitHydrogens(molecule)
        return molecule
    def generateCoordinates(self,molecule):
        sdg = self.cdk.layout.StructureDiagramGenerator(molecule)
        sdg.generateCoordinates()
        return sdg.getMolecule()
    def MolFromSmiles(self,smiles):
        molecule = self.sp.parseSmiles(smiles)
        self.cdk.tools.manipulator.AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule) #
        self.cdk.aromaticity.CDKHueckelAromaticityDetector.detectAromaticity(molecule)
        ha=self.cdk.tools.CDKHydrogenAdder.getInstance(self.builder)
        ha.addImplicitHydrogens(molecule)
        molecule = self.generateCoordinates(molecule)
        return molecule
    def MolToSmiles(self,molecule):
        return self.sg.createSMILES(molecule) #sg.createSMILES(molecule,True)
    def MolToInchiKey(self,molecule):
        igf = self.cdk.inchi.InChIGeneratorFactory.getInstance()
        ig = igf.getInChIGenerator(molecule)
        # print ig
        return ig.getInchiKey()
    def GetExtendedAtomMass(self,atom):
        mass=self.isof.getMajorIsotope(atom.getSymbol()).getExactMass().floatValue()
        try:
            hc=atom.getImplicitHydrogenCount().intValue()
        except:
            hc=0
        return mass+self.Hmass*hc
    def GetFormulaProps(self,mol):
        formula=self.cdk.tools.manipulator.MolecularFormulaManipulator.getMolecularFormula(mol)
        formula_string = self.cdk.tools.manipulator.MolecularFormulaManipulator.getString(formula)
        mim = self.cdk.tools.manipulator.MolecularFormulaManipulator.getMajorIsotopeMass(formula)
        return mim,formula_string
    def FragmentToSmiles(self,mol,atomlist):
        return self.sg.createSMILES(self.acm.extractSubstructure(mol,atomlist))
    def FragmentToInchiKey(self,mol,atomlist):
        ac=self.acm.extractSubstructure(mol,atomlist)
        igf = self.cdk.inchi.InChIGeneratorFactory.getInstance()
        ig = igf.getInChIGenerator(ac)
        return ig.getInchiKey()

Chem = CDKengine()

typew={"AROMATIC":3.0,\
       "DOUBLE":2.0,\
       "TRIPLE":3.0,\
       "SINGLE":1.0}
#typew={Chem.rdchem.BondType.AROMATIC:3.0,\
#         Chem.rdchem.BondType.DOUBLE:2.0,\
#         Chem.rdchem.BondType.TRIPLE:3.0,\
#         Chem.rdchem.BondType.SINGLE:1.0}
ringw={False:1,True:1}
heterow={False:2,True:1}
missingfragmentpenalty=10
max_small_losses=1
weigh_scores=False


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

    def calc_mim(self,mol):
        mim=0
        for atom in range(mol.getAtomCount()):
            mim+=(mims[atom.GetAtomicNum()]+Hmass*\
              (atom.GetNumImplicitHs()+atom.GetNumExplicitHs()))
        return mim

    def add_structure(self,molblock,name,prob,level,sequence,isquery,mim=None,inchikey=None,molform=None,reference=None):
        mol=Chem.MolFromMolBlock(molblock)
        if inchikey == None:
            inchikey=Chem.MolToInchiKey(mol)[:14]
            # inchikey=Chem.MolToSmiles(mol)
        if mim == None or molform == None:
            mim,molform=Chem.GetFormulaProps(mol)
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
            reference=reference
            #logp=Chem.Crippen.MolLogP(m)
            )
        try:
            dupid = self.db_session.query(Metabolite).filter_by(smiles=inchikey).one()
            if dupid.probability < prob:
                self.db_session.delete(dupid)
                metab.metid=dupid.metid
                self.db_session.add(metab)
                sys.stderr.write('Duplicate structure: '+sequence+' '+inchikey+' - old one removed\n')
                # TODO remove any fragments related to this structure as well
            else:
                sys.stderr.write('Duplicate structure: '+sequence+' '+inchikey+' - kept old one\n')
                metab = dupid
        except NoResultFound:
            self.db_session.add(metab)
            sys.stderr.write('Added: '+name+'\n')
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
                    reactionsequence=sequence[:-1]
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

#class ScanType(object):
#    def __init__(self,scanid,mslevel):
#        self.peaks=[]
#        self.scanid=scanid
#        self.mslevel=mslevel
#    
#class PeakType(object):
#    def __init__(self,mz,intensity,scanid,missing_fragment_score):
#        self.mz=mz
#        self.intensity=intensity
#        self.scan=scanid
#        self.childscan=None
#        self.missing_fragment_score=missing_fragment_score
#
##    def massmatch_rel(self,mim,low,high):
##        for x in range(low,high+1):
##            # if self.mz/me.precision < mim+x*Hmass < self.mz*me.precision:
##            if self.mz/1.000005 < mim+x*Hmass < self.mz*1.000005:
##            # if mim/precision < self.mz-x*Hmass < mim*precision:
##                return x
##        else:
##            return False
#
#class HitType(object):
#    def __init__(self,peak,fragment,score,bondbreaks,mass,deltaH):
#        self.mz = peak.mz
#        self.intensity = peak.intensity
#        self.intensity_weight = peak.missing_fragment_score / missingfragmentpenalty
#        self.scan = peak.scan
#        self.fragment = fragment
#        self.score = score
#        self.breaks = bondbreaks
#        self.mass = mass
#        self.deltaH = deltaH
#        self.bonds = []
#        self.allbonds = 0
#        self.besthits=[]
#        self.atomstring=''
#        self.atomlist=[]
#        #print "childscan",peak.childscan

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
            scan.peaks.append(types.PeakType(dbpeak.mz,dbpeak.intensity,scan.scanid,missingfragmentpenalty*(dbpeak.intensity**0.5)))
        dbchildscans=self.db_session.query(Scan).filter(Scan.precursorscanid==scan.scanid).all()
        for dbchildscan in dbchildscans:
            for peak in scan.peaks:
                if -self.precursor_mz_precision < dbchildscan.precursormz-peak.mz < self.precursor_mz_precision:
                    peak.childscan=self.build_spectrum(dbchildscan)
                    for childpeak in peak.childscan.peaks:
                        peak.missing_fragment_score+=childpeak.missing_fragment_score
                    break
            else:
                if dbchildscan.precursorintensity >= cutoff:
                    scan.peaks.append(types.PeakType(dbchildscan.precursormz,dbchildscan.precursorintensity,scan.scanid,missingfragmentpenalty*(dbchildscan.precursorintensity**0.5)))
                    scan.peaks[-1].childscan=self.build_spectrum(dbchildscan)
                    for childpeak in scan.peaks[-1].childscan.peaks:
                        scan.peaks[-1].missing_fragment_score+=childpeak.missing_fragment_score
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
                mass=peak.mz-self.ionisation_mode*Hmass
                if not ((not self.use_all_peaks) and peak.childscan==None):
                    result = c.execute('SELECT * FROM molecules WHERE mim BETWEEN ? AND ?' , (mass/self.precision,mass*self.precision))
                    for (id,mim,molblock,chebi_name) in result:
                        db_candidates[id]=[molblock,chebi_name]
                    print str(mass)+' --> '+str(len(db_candidates))+' candidates'
        return db_candidates

    def get_pubchem_candidates(self):
        dbfilename = '/media/PubChem/Pubchem_MAGMa.db'
        # dbfilename = '/home/ridderl/PRI/test_out.db'
        conn = sqlite3.connect(dbfilename)
        conn.text_factory=str
        c = conn.cursor()
        db_candidates={}
        # First a dictionary is created with candidate molecules based on all level1 peaks
        # In this way duplicates originating from repeated detection of the same component
        # are removed before attempting to add the candidates to the database
        for scan in self.scans:
            for peak in scan.peaks:
                mass=peak.mz-self.ionisation_mode*Hmass
                if not ((not self.use_all_peaks) and peak.childscan==None):
                    result = c.execute('SELECT * FROM molecules WHERE refscore > 2 AND mim BETWEEN ? AND ?' , (mass/self.precision,mass*self.precision))
                    # result = c.execute('SELECT * FROM connectivities JOIN isomers ON isomers.cid = (SELECT cid FROM isomers WHERE isomers.conn_id = connectivities.id LIMIT 1) WHERE connectivities.mim BETWEEN ? AND ?', 
                    #                             (mass/self.precision,mass*self.precision))
                    for (id,cid,mim,molblock,inchikey,molform,name,refscore) in result:
                        db_candidates[id]={'mim':mim,
                                           'mol':zlib.decompress(molblock),
                                           'inchikey':inchikey,
                                           'molform':molform,
                                           'cid':cid,
                                           'name':name,
                                           'refscore':refscore
                                           }
                    print str(scan.scanid)+','+str(peak.mz)+' --> '+str(len(db_candidates))+' candidates'
        return db_candidates

#    def search_all_structures(self):
#        logging.warn('Searching all structures')
#        for structure in self.db_session.query(Metabolite).all():
#            structure.nhits=search_structure(structure,self.scans)

    def search_structures(self,metids=None,ncpus=1):
        logging.warn('Searching structures')
        if metids==None:
            structures = self.db_session.query(Metabolite).all()
        else:
            structures = self.db_session.query(Metabolite).filter(Metabolite.metid.in_(metids)).all()
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
            int_mass=int(round(structure.mim+self.ionisation_mode*Hmass))
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
                   job_server.submit(search_structure,(structure,
                              peaks,
                              self.max_broken_bonds,
                              max_small_losses,
                              self.precision,
                              self.mz_precision_abs,
                              self.use_all_peaks,
                              self.ionisation_mode
                              ),(
                              CDKengine,
                              ),(
                              "numpy",
                              "jpype",
                              "magma.types"
                              )
                           )))
        for structure,job in jobs:
            raw_result=job(raw_result=True)
            hits,sout = pickle.loads(raw_result)
            # print sout
            sys.stderr.write('Metabolite '+str(structure.metid)+': '+str(structure.origin)+'\n')
            structure.nhits=len(hits)
            self.db_session.add(structure)
            for hit in hits:
                sys.stderr.write('Scan: '+str(hit.scan)+' - Mz: '+str(hit.mz)+' - ')
                # storeFragment(metabolite.metid,scan.precursorscanid,scan.precursorpeakmz,2**len(metabolite.atombits)-1,deltaH)
                sys.stderr.write('Score: '+str(hit.score)+'\n')
                # outfile.write("\t"+str(hit.score/fragment_store.get_avg_score()))
                self.store_hit(hit,structure.metid,0)
        self.db_session.commit()

    def store_hit(self,hit,metid,parentfragid):
        global fragid
        fragid+=1
        currentFragid=fragid
        score=hit.score
        deltappm=None
        if score != None:
            score=score/hit.intensity_weight
            deltappm=(hit.mz+hit.deltaH*Hmass-hit.mass)/hit.mz*1e6
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

def search_structure(structure,peaks,max_broken_bonds,max_small_losses,precision,mz_precision_abs,use_all_peaks,ionisation_mode):
    Fragmented=False
    Chem=CDKengine()
    mol=Chem.MolFromMolBlock(str(structure.mol))

    typew={"AROMATIC":3.0,\
           "DOUBLE":2.0,\
           "TRIPLE":3.0,\
           "SINGLE":1.0}
    #typew={Chem.rdchem.BondType.AROMATIC:3.0,\
    #         Chem.rdchem.BondType.DOUBLE:2.0,\
    #         Chem.rdchem.BondType.TRIPLE:3.0,\
    #         Chem.rdchem.BondType.SINGLE:1.0}
    global missingfragmentpenalty
    ringw={False:1,True:1}
    heterow={False:2,True:1}
    missingfragmentpenalty=10
    max_small_losses=1
    weigh_scores=False
    
    
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

    class GrowingEngine(object):
        def __init__(self,mol):
#                self.natoms=mol.GetNumAtoms() # RDKit
            self.natoms=mol.getAtomCount()
            self.FragmentBreaks={}
            self.FragmentMass={}
            self.atombits=[]          # [1,2,4,8,16,....]
            self.atommass={}          # {atombit:atommass(incl. attached hydrogens)}
            self.bondbits=[]
            self.bondscore=[]

            for x in range(self.natoms):
                self.atombits.append(2**x)
#                    an=mol.GetAtomWithIdx(x).GetAtomicNum() #RDKit
                an=mol.getAtom(x).getAtomicNumber()
                self.atommass[2**x]=Chem.GetExtendedAtomMass(mol.getAtom(x))
#                for x in mol.GetBonds():
            for x in range(mol.getBondCount()):
                bondbit=0
                for a in mol.getBond(x).atoms().iterator():
                    bondbit=bondbit | mol.getAtomNumber(a)
                self.bondbits.append(bondbit)
                
#                    self.bondbits.append(self.atombits[x.GetBeginAtomIdx()]|\
#                                         self.atombits[x.GetEndAtomIdx()])
#                    self.bondscore.append(typew[x.GetBondType()]*ringw[x.IsInRing()]*\
#                                          heterow[x.GetBeginAtom().GetAtomicNum() != 6 or \
#                                                     x.GetEndAtom().GetAtomicNum() != 6])

        def grow(self,fragment):
            NewFragments=[]
            for bond in self.bondbits:
                if 0 < (fragment & bond) < bond:
                    NewFragments.append(fragment|bond)
            self.FragmentBreaks[fragment]=len(NewFragments)
            if len(NewFragments) <= max_broken_bonds+3:
                for NewFragment in NewFragments:
                    if NewFragment not in self.FragmentMass:
                        self.FragmentMass[NewFragment]=self.FragmentMass[fragment]+self.atommass[NewFragment^fragment]
                        self.grow(NewFragment)

        def generate_fragments(self):
            for atom in self.atombits:
                self.FragmentMass[atom]=self.atommass[atom]
                self.grow(atom)
            sys.stderr.write('N fragments created: '+str(len(self.FragmentMass))+"\n")

            # first items fragment_masses and fragment_info represent the complete molecule
            # represent the complete molecule
            fragment_store=FragmentStore()
            for fragment in self.FragmentMass:
                if self.FragmentBreaks[fragment]<=max_broken_bonds:
                    score=0
                    for bond in range(len(self.bondbits)):
                        if 0 < (fragment & self.bondbits[bond]) < self.bondbits[bond]:
                            score += self.bondscore[bond]
                    fragment_store.add_fragment(fragment,self.FragmentMass[fragment],score,self.FragmentBreaks[fragment])
            return fragment_store


    class FragmentationEngine(object):
        def __init__(self,mol):
            self.natoms=mol.getAtomCount()  # number of atoms in the molecule
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
                self.atom_masses.append(Chem.GetExtendedAtomMass(mol.getAtom(x)))
                if mol.getAtom(x).getSymbol() == 'O' and mol.getAtom(x).getImplicitHydrogenCount() == 1:
                    self.neutral_loss_atoms.append(x)
            for x in range(mol.getBondCount()):
                bond=0
                heterobond=False
                for a in mol.getBond(x).atoms().iterator():
                    self.bonded_atoms[mol.getAtomNumber(a)].append(mol.getAtomNumber(mol.getBond(x).getConnectedAtom(a)))
                    bond = bond | 1<<mol.getAtomNumber(a)
                    if a.getSymbol() != 'C':
                        heterobond=True
                if mol.getBond(x).getFlag(4) == 1:
                    bondscore = typew["AROMATIC"]
                else:
                    bondscore = typew[mol.getBond(x).getOrder().toString()]
                bondscore*=heterow[heterobond]
                    

#                    for bond in mol.GetAtomWithIdx(x).GetBonds():
#                        self.bonded_atoms[-1].append(bond.GetBeginAtomIdx()+bond.GetEndAtomIdx()-x)
#                    atom = mol.GetAtomWithIdx(x)
#                    self.atom_masses.append(mims[atom.GetAtomicNum()]+Hmass*(atom.GetNumImplicitHs()+atom.GetNumExplicitHs()))
#                    if atom.GetAtomicNum()==8 and (atom.GetNumImplicitHs()+atom.GetNumExplicitHs())==1:
#                        self.neutral_loss_atoms.append(x)
#                for x in mol.GetBonds():
#                    bond = (1<<x.GetBeginAtomIdx()) | (1<<x.GetEndAtomIdx())
#                    bondscore = typew[x.GetBondType()]*ringw[x.IsInRing()]*\
#                                          heterow[x.GetBeginAtom().GetAtomicNum() != 6 or \
#                                                     x.GetEndAtom().GetAtomicNum() != 6]
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
            for step in range(max_broken_bonds):                    # perform fragmentation for nstep steps
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
                                    if bondbreaks<=max_broken_bonds and score < (missingfragmentpenalty+5):
                                        self.new_fragments.add(frag)
                                        self.total_fragments.add(frag)
                                        self.add_fragment(frag,self.calc_fragment_mass(frag),score,bondbreaks)
                self.current_fragments=self.new_fragments
                self.new_fragments=set([])
            for step in range(max_small_losses):                    # number of OH losses
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
                           numpy.hstack((numpy.zeros(max_broken_bonds+max_small_losses-bondbreaks),
                                      numpy.arange(-bondbreaks-1,bondbreaks+2)*Hmass+fragmentmass,
                                      numpy.zeros(max_broken_bonds+max_small_losses-bondbreaks)
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

        def find_hit(self,peak,parent):
            # subfrag=(1<<mol.GetNumAtoms())-1 # remove hierarchical constraint
            # besthit=None
            besthit=gethit(peak,0,None,0,0,0)
            # print peak.missingfragmentscore
            # result=numpy.where(numpy.where(self.fragment_masses < (peak.mz+mz_precision),self.fragment_masses,0) > (peak.mz-mz_precision))
            result=numpy.where(numpy.where(self.fragment_masses < max(peak.mz*precision,peak.mz+mz_precision_abs),
                                     self.fragment_masses,0) > min(peak.mz/precision,peak.mz-mz_precision_abs))
            for i in range(len(result[0])):
                fid=result[0][i]
                if self.fragments[fid] & parent == self.fragments[fid]:
                # if self.info[fid][0] & subfrag == self.info[fid][0]:
                    hitscore=self.scores[fid]*(peak.intensity**0.5) # *(peak.mz**3) # <--
                    # hitscore=self.score_fragment(self.fragments[fid],subfrag) # *(peak.intensity**0.5) # <--
                    # losses=numpy.nonzero(self.fragments==subfrag^self.fragments[fid])
                    # if len(losses[0])==1:
                    #     hitscore=self.scores[losses[0][0]]
                    hit=gethit(peak,self.fragments[fid],hitscore,self.bondbreaks[fid],
                                self.fragment_masses[fid][max_broken_bonds+max_small_losses+1],max_broken_bonds+max_small_losses+1-result[1][i])
                    if besthit.score==None or besthit.score > hit.score or \
                           (besthit.score == hit.score and abs(besthit.deltaH) > abs(hit.deltaH)) or \
                           self.score_fragment_rel2parent(besthit.fragment,parent) > self.score_fragment_rel2parent(hit.fragment,parent):
                        besthit=hit
            return besthit

#            def score_fragment(self,fragment,subfrag):
#                score=0.0
#                for bond in self.bonds:
#                    if 0 < (fragment & bond) < (subfrag & bond):
#                    #if 0 < (fragment & bond) < bond:
#                        score+=self.bondscore[bond]
#                return score


    def massmatch(peak,mim,low,high):
        for x in range(low,high+1):
            #if self.mz-me.mz_precision < mim+x*Hmass < self.mz+me.mz_precision:
            if peak.mz/precision < mim+x*Hmass < peak.mz*precision:
                return x
        else:
            return False

    #def findhit(self,childscan,parent):
    def gethit (peak,fragment,score,bondbreaks,mass,deltaH):
        try:
            hit=types.HitType(peak,fragment,score,bondbreaks,mass,deltaH)
        except:
            hit=magma.types.HitType(peak,fragment,score,bondbreaks,mass,deltaH)
        #hit.__module__="magma"
        #rint hit.__module__
        #print os.path.splitext(os.path.basename(__file__))[0]
        if fragment>0 and peak.childscan!=None and len(peak.childscan.peaks) > 0: # fragment=0 means it is a missing fragment
            n_child_peaks=len(peak.childscan.peaks)
            # self.score = self.score/(n_child_peaks+1) + self.findhits(peak.childscan,self.fragment)*n_child_peaks/(n_child_peaks+1)
            total_score=0.0
            total_count=0.0
            for childpeak in peak.childscan.peaks:
                if weigh_scores:
                    weight=(childpeak.intensity**0.5) # *(peak.mz**3)
                else:
                    weight=1
                total_count+=weight
                besthit=fragment_engine.find_hit(childpeak,fragment)
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
        atomlist=[]
        atomstring=''
        if hit.fragment != 0:
            for atom in range(mol.getAtomCount()):
                if ((1<<atom) & hit.fragment):
                    atomstring+=','+str(atom)
                    atomlist.append(atom)
            hit.atomstring=atomstring
            hit.atomlist=atomlist
            try:
                hit.inchikey=Chem.FragmentToInchiKey(mol,hit.atomlist)[:14]
            except:
                exit('failted inchi for: '+atomstring+'--'+str(hit.fragment))
            if len(hit.besthits)>0:
                for childhit in hit.besthits:
                    if childhit != None: # still need to work out how to deal with missed fragments
                        add_fragment_data_to_hit(childhit)

# for peak in self.db_session.query(Peak).filter(Scan.mslevel)==1.filter(Peak.intensity>MSfilter)
    hits=[]
    for peak in peaks:
        if not ((not use_all_peaks) and peak.childscan==None):
            protonation=ionisation_mode-(structure.molformula.find('+')>=0)*1
            deltaH=massmatch(peak,structure.mim,protonation,protonation)
            if type(deltaH)==int:
                if not Fragmented:
                    #sys.stderr.write('\nMetabolite '+str(structure.metid)+': '+str(structure.origin)+' '+str(structure.reactionsequence)+'\n')
                    #sys.stderr.write('Mim: '+str(structure.mim)+'\n')
                    fragment_engine=FragmentationEngine(mol)
                    #fragment_engine=GrowingEngine(mol)
                    fragment_engine.generate_fragments()
                    #sys.stderr.write('N fragments kept: '+str(len(fragment_engine.fragments))+"\n")
                    Fragmented=True
                hit=gethit(peak,(1<<mol.getAtomCount())-1,0,0,structure.mim,-deltaH)
                add_fragment_data_to_hit(hit)
                hits.append(hit)
    return hits



