#!/usr/bin/env python

import sys
import base64
import time
import re
import os
import sqlite3
import struct
import zlib
import copy
import pkg_resources
import logging
import json
from lxml import etree
from sqlalchemy import create_engine, desc
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import func
from models import Base, Molecule, Reaction, fill_molecules_reactions, Scan, Peak, Fragment, Run
import requests
import macauthlib  # required to update callback url
from requests.auth import AuthBase
import pp
import cPickle as pickle
import types
from magma.errors import FileFormatError, DataProcessingError
import pars

import ConfigParser
config = ConfigParser.ConfigParser()
# read config file from current working directory or users home dir
config.read(['magma_job.ini', os.path.expanduser('~/magma_job.ini')])

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms, Descriptors
from rdkit.rdBase import DisableLog
DisableLog('rdApp.warning')

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger('MagmaLogger')


class MagmaSession(object):

    """ MAGMa session object, which is linked to a sqlite3 db_name.db file
        and enables starting different engines to perform MAGMa functions """

    def __init__(self, db_name, description="", loglevel='info'):
        if db_name == None:
            engine = create_engine('sqlite://', echo=False)
        else:
            engine = create_engine('sqlite:///' + db_name, echo=False)
        session = sessionmaker()
        session.configure(bind=engine)
        self.db_session = session()
        Base.metadata.create_all(engine)
        try:
            rundata = self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.description is None:
            rundata.description = unicode(description)
        self.db_session.add(rundata)
        self.db_session.commit()
        logger.setLevel(getattr(logging, loglevel.upper()))

    def get_structure_engine(self, pubchem_names=False, call_back_url=None):
        return StructureEngine(self.db_session, pubchem_names, call_back_url)

    def get_ms_data_engine(self,
                           ionisation_mode=1,
                           abs_peak_cutoff=1000,
                           mz_precision=5.0,
                           mz_precision_abs=0.001,
                           precursor_mz_precision=0.005,
                           max_ms_level=10,
                           call_back_url=None
                           ):
        return MsDataEngine(self.db_session, ionisation_mode, abs_peak_cutoff, mz_precision, mz_precision_abs,
                            precursor_mz_precision, max_ms_level, call_back_url)

    def get_annotate_engine(self,
                            skip_fragmentation=False,
                            max_broken_bonds=3,
                            max_water_losses=1,
                            ms_intensity_cutoff=1e6,
                            msms_intensity_cutoff=5,
                            use_all_peaks=False,
                            adducts=None,
                            max_charge=1,
                            call_back_url=None
                            ):
        return AnnotateEngine(self.db_session, skip_fragmentation, max_broken_bonds, max_water_losses,
                              ms_intensity_cutoff, msms_intensity_cutoff, use_all_peaks, adducts, max_charge, call_back_url)

    def get_export_molecules_engine(self):
        return ExportMoleculesEngine(self.db_session)

    def get_call_back_engine(self, id, key):
        return CallBackEngine(id, key)

    def fill_molecules_reactions(self):
        fill_molecules_reactions(self.db_session)

    def commit(self):
        self.db_session.commit()

    def close(self):
        self.db_session.close()


class CallBackEngine(object):

    """ Callback engine object, used to send job status info to callback url"""

    def __init__(self, url):
        self.access_token = config.get('magma job', 'macs.id')
        self.mac_key = config.get('magma job', 'macs.key')
        self.url = url
        # send update to call_back_url every 2 seconds
        self.update_interval = 2
        self.update_time = time.time()

    def update_callback_url(self, status, elapsed_time=None, time_limit=None, force=False):
        class HTTPMacAuth(AuthBase):

            """Attaches HTTP Basic Authentication to the given Request object."""

            def __init__(self, id, key):
                self.id = id
                self.key = key

            def __call__(self, r):
                r.headers['Authorization'] = macauthlib.sign_request(
                    r, id=self.id, key=self.key)
                return r

        update_last = time.time() - self.update_time
        # update status every second
        if force or update_last > self.update_interval:
            if elapsed_time is not None:
                status += '<h3>Time: %02d:%02d:%02d' % (
                    elapsed_time // 3600, (elapsed_time % 3600) // 60, elapsed_time % 60)
                if time_limit is not None:
                    status += ' / max. %02d:%02d:00 (%d%%)' % (
                        time_limit // 60, time_limit % 60, elapsed_time / time_limit / 60 * 100)
                status += '</h3>'
            self.update_time = self.update_time + \
                update_last // self.update_interval * self.update_interval
            requests.put(
                self.url, status, auth=HTTPMacAuth(self.access_token, self.mac_key))


class MetabolizeEngine(object):

    """ Engine to perform in silico reactions on chemical structures """

    def __init__(self):
        metabolism_files = {
            "phase1": pkg_resources.resource_filename(  # @UndefinedVariable
                'magma', "data/phase1.smirks"),
            "phase2": pkg_resources.resource_filename(  # @UndefinedVariable
                'magma', "data/phase2.smirks"),
            "phase1_selected": pkg_resources.resource_filename(  # @UndefinedVariable
                'magma', "data/phase1_selected.smirks"),
            "phase2_selected": pkg_resources.resource_filename(  # @UndefinedVariable
                'magma', "data/phase2_selected.smirks"),
            "gut": pkg_resources.resource_filename(  # @UndefinedVariable
                'magma', "data/gut.smirks"),
            "glycosidase": pkg_resources.resource_filename(  # @UndefinedVariable
                'magma', "data/glycosidase.smirks"),
        }
        self.reactions = {}
        for m in metabolism_files:
            self.reactions[m] = []
            f = open(metabolism_files[m])
            for line in f:
                if line[0] != "#" and len(line) > 1:
                    splitline = line[:-1].split()
                    self.reactions[m].append(
                        (splitline[1], AllChem.ReactionFromSmarts(splitline[0])))

    def react(self, reactant, reaction):
        """ Apply reaction to reactant and return products """
        ps = reaction.RunReactants([reactant])
        products = []
        for product in ps:
            frags = (Chem.GetMolFrags(product[0], asMols=True))
            for p in frags:
                q = copy.copy(p)
                Chem.SanitizeMol(q)
                self.gen_coords(q)
                products.append(q)
        return products

    def gen_coords(self, mol):
        """ Calculate 2D positions for atoms without coordinates """
        conf = mol.GetConformer(0)
        coordDict = {}
        # Put known coordinates in coordDict
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            if pos.x != 0.0 or pos.y != 0.0:
                coordDict[i] = Geometry.Point2D(pos.x, pos.y)
        if len(coordDict) > 0:
            # calculate average length of all bonds with coordinates
            total = 0
            n = 0
            for bond in mol.GetBonds():
                b = bond.GetBeginAtomIdx()
                e = bond.GetEndAtomIdx()
                if b in coordDict and e in coordDict:
                    n += 1
                    total += rdMolTransforms.GetBondLength(conf, b, e)
            av = total / n
            # compute coordinates for new atoms, keeping known coordinates
            AllChem.Compute2DCoords(mol, coordMap=coordDict, bondLength=av)

    def metabolize(self, mol, metabolism, endpoints=False):
        """ Metabolize a molecule based on a set of metabolism rules and return products"""
        metabolites = {}  # {ikey:[rulename,product mol]}
        if metabolism not in self.reactions:
            return metabolites
        # if endpoints, name refers to the metabolism type
        name = metabolism
        for reaction in self.reactions[metabolism]:
            if not endpoints:
                # else name refers to the rule generating a product
                name = reaction[0]
            products = self.react(mol, reaction[1])
            for p in products:
                ikey = AllChem.InchiToInchiKey(AllChem.MolToInchi(p))[:14]
                if ikey not in metabolites:
                    if endpoints:
                        # Try further reactions, otherwise store end product
                        endproducts = self.metabolize(p, metabolism, True)
                        if len(endproducts) > 0:
                            metabolites.update(endproducts)
                        else:
                            metabolites[ikey] = [name, p]
                    else:
                        metabolites[ikey] = [name, p]
        return metabolites


def get_molecule(molblock, name, refscore, predicted, mim=None, smiles=None, natoms=None,
                 inchikey14=None, molform=None, reference=None, logp=None, mass_filter=9999):
    """ Returns a Molecule with the given attributes """
    if inchikey14 is None or mim is None or molform is None or logp is None or natoms is None:
        mol = Chem.MolFromMolBlock(molblock)
        if mol is None:
            return
        inchikey14 = Chem.AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))[:14]
        smiles = Chem.MolToSmiles(mol)
        # Calc mim
        mim = 0.0
        for atom in mol.GetAtoms():
            mim += pars.mims[atom.GetSymbol()] + pars.Hmass * \
                (atom.GetNumImplicitHs() + atom.GetNumExplicitHs())
        molform = Chem.rdMolDescriptors.CalcMolFormula(mol)
        natoms = mol.GetNumHeavyAtoms()
        logp = Chem.Crippen.MolLogP(mol)
    if mim > mass_filter:
        return
    return Molecule(
        mol=unicode(molblock),
        refscore=refscore,
        inchikey14=unicode(inchikey14),
        smiles=unicode(smiles),
        formula=unicode(molform),
        predicted=predicted,
        name=unicode(name),
        nhits=0,
        mim=mim,
        natoms=natoms,
        reference=unicode(reference),
        logp=logp,
        reactionsequence={})


class StructureEngine(object):

    """Engine to add and query MAGMa candidate molecules"""

    def __init__(self, db_session, pubchem_names=False, call_back_url=None):
        self.db_session = db_session
        self.pubchem_names = pubchem_names
        if pubchem_names:
            self.pubchem_engine = PubChemEngine('pubchem')
        self.metabolize_engine = None
        if call_back_url is not None:
            self.call_back_engine = CallBackEngine(call_back_url)
        else:
            self.call_back_engine = None

    def read_sdf(self, file_name, mass_filter):
        logger.info('READING SDF (Structure Data File)')
        molids = set([])
        for mol in Chem.SDMolSupplier(file_name):
            molids.add(self.add_structure(Chem.MolToMolBlock(mol), mol.GetProp(
                '_Name'), None, predicted=0, mass_filter=mass_filter))
        logger.info(str(len(molids)) + ' molecules added to library\n')

    def read_smiles(self, smiles, mass_filter=9999):
        """ Read smiles file,
            each line consists of a smiles string
            optionally followed by a name in which spaces a replaced by underscores """
        logger.info('READING SMILES')
        try:
            smiles_file = open(smiles)
        except:
            smiles_file = [smiles]
        molids = set([])
        for line in smiles_file:
            line = line.strip()
            if line == "":
                continue
            splitline = line.split()
            smiles = splitline[0]
            if len(splitline) > 1:
                name = '_'.join(splitline[1:])
            else:
                name = ''
            try:
                mol = Chem.MolFromSmiles(smiles)
                mol.SetProp('_Name', name)
                AllChem.Compute2DCoords(mol)
                molids.add(self.add_structure(
                    Chem.MolToMolBlock(mol), name, None, predicted=0, mass_filter=mass_filter))
            except:
                logger.warn(
                    'Failed to read smiles: ' + smiles + ' (' + name + ')')
        logger.info(str(len(molids)) + ' molecules added to library\n')

    def add_structure(self, molblock, name, refscore, predicted,
                      mim=None, smiles=None, natoms=None, inchi=None, molform=None, reference=None, logp=None, mass_filter=9999):
        """ Get molecule with given attributes and add to database
            Return molid """
        molecule = get_molecule(
            molblock, name, refscore, predicted, mim, natoms, inchi, molform, reference, logp)
        return self.add_molecule(molecule, mass_filter)

    def add_molecule(self, newmol, mass_filter=9999, check_duplicates=True, merge=False):
        """ Add molecule to database
            Optionally lookup name in pubchem
            Return molid """
        # Checking duplicates is slow, option to turn off should only be used
        # when candidates from a single chemical db query are added
        if check_duplicates:
            dup = self.db_session.query(Molecule).filter_by(
                inchikey14=newmol.inchikey14).first()
            if dup is not None:
                if merge:
                    if dup.name == "" and newmol.name != "":
                        dup.name = newmol.name
                    if dup.reference == "None" and newmol.reference != "":
                        dup.reference = newmol.reference
                    if newmol.refscore > 0 and newmol.refscore > dup.refscore:
                        dup.refscore = newmol.refscore
                    newmol = dup
                else:
                    if dup.refscore < newmol.refscore:
                        newmol.molid = dup.molid
                        self.db_session.delete(dup)
                        logger.info('Duplicate structure, first one removed: ' + newmol.name)
                    # TODO remove any fragments related to this structure as well
                    else:
                        logger.info(
                            'Duplicate structure, kept new one: ' + newmol.name)
                        return
        if self.pubchem_names:
            in_pubchem = self.pubchem_engine.check_inchi(newmol.mim, newmol.inchikey14)
            if in_pubchem != False:
                name, reference, refscore = in_pubchem
                if newmol.name == '':
                    newmol.name = unicode(name, 'utf-8', 'xmlcharrefreplace')
                if newmol.reference == "None":
                    newmol.reference = unicode(reference)
                if newmol.refscore is None:
                    newmol.refscore = refscore
        self.db_session.add(newmol)
        logger.debug('Added molecule: ' + newmol.name)
        self.db_session.flush()
        return newmol.molid

    def metabolize(self, molid, metabolism, endpoints=False):
        """ Metabolize a molecule based on a set of metabolism rules and store
            reactions and products in reactions and molecules tables
            Return molids """
        try:
            p = self.db_session.query(Molecule).filter_by(molid=molid).one()
        except:
            logger.warn('Molecule record %s does not exist.', molid)
            return
        parent = Chem.MolFromMolBlock(p.mol)
        if self.metabolize_engine is None:
            self.metabolize_engine = MetabolizeEngine()
        molids = set()
        products = self.metabolize_engine.metabolize(parent, metabolism, endpoints)
        for product in products.values():
            molecule = get_molecule(Chem.MolToMolBlock(product[1]), "", None, 1)
            new_molid = self.add_molecule(molecule, merge=True)
            molids.add(new_molid)
            reactid = self.db_session.query(Reaction.reactid).filter(
                            Reaction.reactant == molid,
                            Reaction.product == new_molid,
                            Reaction.name == unicode(product[0])).all()
            if len(reactid) == 0:
                react = Reaction(
                    reactant=molid,
                    product=new_molid,
                    name=unicode(product[0])
                )
                self.db_session.add(react)
                self.db_session.flush()
        if endpoints and len(molids) == 0:  # molid itself is endpoint
            molids.add(molid)
        return molids

    def metabolize_all(self, metabolism, endpoints=False):
        """ Metabolize all current molecules based on a set of metabolism rules """
        logger.info('Metabolize all')
        parentids = self.db_session.query(Molecule.molid).all()
        molids = set([])
        for parentid, in parentids:
            molids |= self.metabolize(parentid, metabolism, endpoints)
        return molids

    def run_scenario(self, scenario, time_limit=None):
        """ Metabolize according to scenario """
        logger.info('RUNNING METABOLIC SCENARIO')
        if time_limit is None:
            result = self.db_session.query(Molecule.molid).all()
            molids = {x[0] for x in result}
        else:
            # in case of time_limit only metabolize the first compound
            molids = set([self.db_session.query(Molecule.molid).first()[0]])
        start_time = time.time()
        for step in range(len(scenario)):
            action, value = scenario[step]
            logger.info("Stage " + str(step + 1) + ":")
            endpoints = False
            prev_molids = molids
            if action == 'mass_filter':
                result = self.db_session.query(Molecule.molid).filter(
                    Molecule.mim < float(value), Molecule.molid.in_(molids)).all()
                molids = {x[0] for x in result}
                logger.info('    from ' + str(len(prev_molids)) + ' compounds, ' +
                            str(len(molids)) + ' were selected with mass <' + str(value))
            else:
                if value == 'complete':
                    endpoints = True
                    value = 1
                new_molids = set()
                for molid in molids:
                    new_molids |= self.metabolize(molid, action, endpoints)
                    elapsed_time = time.time() - start_time
                    if self.call_back_engine is not None:
                        status = 'Transformation: %s, step 1<br>Metabolites generated: %d' % (
                            action, len(prev_molids) + len(molids) + len(new_molids))
                        self.call_back_engine.update_callback_url(status, elapsed_time, time_limit)
                    if time_limit and elapsed_time > time_limit * 60:
                        break
                else:
                    active_molids = new_molids.difference(molids)
                    molids = new_molids
                    for i in range(1, int(value)):
                        new_molids = set()
                        for molid in active_molids:
                            new_molids |= self.metabolize(molid, action, endpoints)
                            elapsed_time = time.time() - start_time
                            if self.call_back_engine is not None:
                                status = 'Transformation: %s, step %d<br>Metabolites generated: %d' % (
                                    action, i + 1, len(prev_molids) + len(molids) + len(new_molids))
                                self.call_back_engine.update_callback_url(status, elapsed_time, time_limit)
                            if time_limit and elapsed_time > time_limit * 60:
                                break
                        active_molids = new_molids.difference(molids)
                        molids |= new_molids
                        if time_limit and time.time() - start_time > time_limit * 60:
                            break
                logger.info('    ' + str(len(prev_molids)) + ' compounds were metabolized according to ' + \
                            str(action) + ' rules (' + str(int(value)) + ' steps)')
                logger.info('    yielding ' + str(len(molids)) + ' metabolites')
                if not endpoints:
                    molids |= prev_molids
            if time_limit and time.time() - start_time > time_limit * 60:
                if self.call_back_engine is not None:
                    self.call_back_engine.update_callback_url(
                        'Transformation stopped: time limit exceeded', force=True)
                logger.warn('Transformation stopped: time limit exceeded')
                break
        else:
            if self.call_back_engine is not None:
                self.call_back_engine.update_callback_url('Transformations completed', force=True)
        logger.info(str(self.db_session.query(Molecule).count()) + ' molecules in library\n')


class MsDataEngine(object):

    """ Engine to read MS/MS data """

    def __init__(self, db_session, ionisation_mode, abs_peak_cutoff, mz_precision, mz_precision_abs,
                 precursor_mz_precision, max_ms_level, call_back_url=None):
        self.db_session = db_session
        # a small mz_precision_abs is required,
        # even when matching theoretical masses,
        # because of finite floating point precision
        mz_precision_abs = max(mz_precision_abs, 0.000001)
        precursor_mz_precision = max(precursor_mz_precision, 0.000001)
        try:
            rundata = self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.ionisation_mode is None:
            rundata.ionisation_mode = ionisation_mode
        if rundata.abs_peak_cutoff is None:
            rundata.abs_peak_cutoff = abs_peak_cutoff
        if rundata.max_ms_level is None:
            rundata.max_ms_level = max_ms_level
        if rundata.mz_precision is None:
            rundata.mz_precision = mz_precision
        if rundata.mz_precision_abs is None:
            rundata.mz_precision_abs = mz_precision_abs
        if rundata.precursor_mz_precision is None:
            rundata.precursor_mz_precision = precursor_mz_precision
        self.db_session.add(rundata)
        self.db_session.commit()
        if rundata.ionisation_mode == -1:
            self.polarity = "-"
        else:
            self.polarity = "+"
        self.abs_peak_cutoff = rundata.abs_peak_cutoff
        self.max_ms_level = rundata.max_ms_level
        # Difference between lowest and highest possible value
        self.mz_precision = rundata.mz_precision * 2
        # Difference between lowest and highest possible value
        self.precision = 1 + 2 * rundata.mz_precision / 1e6
        self.mz_precision_abs = rundata.mz_precision_abs * 2
        self.precursor_mz_precision = rundata.precursor_mz_precision
        if call_back_url is not None:
            self.call_back_engine = CallBackEngine(call_back_url)
        else:
            self.call_back_engine = None

    def store_mzxml_file(self, mzxml_file, scan_filter=None, time_limit=None):
        """ Read mzxml_file and store scans and peaks in the scans and peaks tables """
        logger.info('READING MZXML FILE')
        rundata = self.db_session.query(Run).one()
        if rundata.ms_filename is not None:
            raise DataProcessingError('Attempt to read MS data twice')
        rundata.ms_filename = unicode(mzxml_file)
        self.db_session.add(rundata)
        self.ms_filename = mzxml_file
        tree = etree.parse(mzxml_file)
        root = tree.getroot()
        namespace = '{' + root.nsmap[None] + '}'
        mzxml_query = namespace + "msRun/" + namespace + "scan"
        prec_scans = [] # in case of non-hierarchical mzXML, also find child scans of scan_filter
        start_time = time.time()
        for mzxmlScan in root.findall(mzxml_query):
            elapsed_time = time.time() - start_time
            if mzxmlScan.attrib['polarity'] == self.polarity and \
                    (scan_filter is None or \
                     mzxmlScan.attrib['num'] == scan_filter or \
                     (int(mzxmlScan.attrib['msLevel'])>1 and mzxmlScan.find(namespace+'precursorMz').attrib['precursorScanNum'] in prec_scans)
                     ):
                self.store_mzxml_scan(mzxmlScan, 0, namespace)
                prec_scans.append(mzxmlScan.attrib['num']) # in case of non-hierarchical mzXML, also find child scans
                if self.call_back_engine is not None:
                    status = 'Reading mzXML, scan: %s' % (prec_scans[-1],)
                    self.call_back_engine.update_callback_url(status, elapsed_time, time_limit)
            if time_limit and elapsed_time > time_limit * 60:
                if self.call_back_engine is not None:
                    self.call_back_engine.update_callback_url('Reading mzXML stopped: time limit exceeded', force=True)
                logger.warn('Reading mzXML stopped: time limit exceeded\n')
                break
        else:
            if self.call_back_engine is not None:
                self.call_back_engine.update_callback_url('Reading mzXML completed', force=True)
        self.db_session.commit()
        logger.info(str(self.db_session.query(Scan).count()) + ' spectra read from file\n')

    def store_mzxml_scan(self, mzxmlScan, precScan, namespace):
        if mzxmlScan.attrib['peaksCount'] == '0':
            return
        try:
            lowmz = float(mzxmlScan.attrib['lowMz'])
            highmz = float(mzxmlScan.attrib['highMz'])
        except:
            lowmz = None
            highmz = None
        try:
            totioncurrent = float(mzxmlScan.attrib['totIonCurrent'])
        except:
            totioncurrent = None
        scan = Scan(
            scanid=int(mzxmlScan.attrib['num']),
            mslevel=int(mzxmlScan.attrib['msLevel']),
            rt=float(mzxmlScan.attrib['retentionTime'].strip('PTS')) / 60,
            lowmz=lowmz,
            highmz=highmz,
            basepeakmz=float(mzxmlScan.attrib['basePeakMz']),
            basepeakintensity=float(mzxmlScan.attrib['basePeakIntensity']),
            totioncurrent=totioncurrent,
            precursorscanid=precScan
        )
        logger.debug('Processing scan ' + str(scan.scanid) + ' (level ' + str(scan.mslevel) + ')')
        comp = None
        for child in mzxmlScan:
            if child.tag == namespace + 'precursorMz':
                scan.precursormz = float(child.text)
                scan.precursorintensity = float(child.attrib['precursorIntensity'])
                if child.attrib.get('precursorScanNum') is not None:
                    # Non-hierarchical mzXML!
                    scan.precursorscanid = float(child.attrib['precursorScanNum'])
                if scan.precursorscanid == 0:
                    # MS>1 scan is not in a hierarchy and does not contain precursor scan information
                    # assume that the preceding (based on retention time) scan at the precursor MS level is the precursor scan
                    scan.precursorscanid = self.db_session.query(Scan.scanid).filter(
                        Scan.mslevel == scan.mslevel - 1).filter(Scan.rt < scan.rt).order_by(desc(Scan.rt)).first()[0]
                    logger.info('Assigning precursor scanid ' + str(scan.precursorscanid) + ' to scan ' + str(scan.scanid) + '\n')
                # check for existing scan with the same precursor
                comp = self.db_session.query(Scan).filter(Scan.precursorscanid == scan.precursorscanid,
                                                          Scan.precursormz == scan.precursormz,
                                                          Scan.precursorintensity == scan.precursorintensity).all()
            if child.tag == namespace + 'peaks':
                decoded = base64.decodestring(child.text)
                try:
                    if child.attrib['compressionType'] == 'zlib':
                        decoded = zlib.decompress(decoded)
                except:
                    pass
                if comp is None or len(comp) == 0:
                    self.store_mzxml_peaks(scan, decoded)
                else:
                    # generate composite spectrum with the existing scan of the same precursor
                    self.merge_spectrum(comp[0], scan, decoded)
            if child.tag == namespace + 'scan' and int(child.attrib['msLevel']) <= self.max_ms_level:
                self.store_mzxml_scan(child, scan.scanid, namespace)
        if comp is None or len(comp) == 0:
            self.db_session.add(scan)
        self.db_session.flush()

    def merge_spectrum(self, existing_scan, newscan, decoded):
        """ Generate composite spectrum: in case of matching m/z values, keep the peak with
            highest intensity """
        logger.info('Merging scans ' + str(existing_scan.scanid) + ' and ' + str(newscan.scanid))
        if existing_scan.lowmz > newscan.lowmz:
            existing_scan.lowmz = newscan.lowmz
        if existing_scan.highmz < newscan.highmz:
            existing_scan.highmz = newscan.highmz
        if existing_scan.basepeakintensity < newscan.basepeakintensity:
            existing_scan.basepeakintensity = newscan.basepeakintensity
            existing_scan.basepeakmz = newscan.basepeakmz
        self.db_session.add(existing_scan)
        tmp_size = len(decoded) / 4
        unpack_format1 = ">%df" % tmp_size
        unpacked = struct.unpack(unpack_format1, decoded)
        for mz, intensity in zip(unpacked[::2], unpacked[1::2]):
            if intensity > self.abs_peak_cutoff:
                matching_peaks = self.db_session.query(Peak).filter(Peak.scanid == existing_scan.scanid,
                                       Peak.mz.between(min(mz / self.precision, mz - self.mz_precision_abs),
                                                       max(mz * self.precision, mz + self.mz_precision_abs)
                                                       )).all()
                if len(matching_peaks) == 0:
                    self.db_session.add(Peak(scanid=existing_scan.scanid, mz=mz, intensity=intensity))
                else:
                    replace = True
                    # Compare intensity of peak with all matching peaks. Of all those,
                    # keep only the one with highest intensity
                    for p in matching_peaks:
                        if intensity > p.intensity:
                            self.db_session.delete(p)
                        else:
                            replace = False
                            intensity = p.intensity
                    if replace:
                        self.db_session.add(Peak(scanid=existing_scan.scanid, mz=mz, intensity=intensity))

    def store_mzxml_peaks(self, scan, decoded):
        tmp_size = len(decoded) / 4
        unpack_format1 = ">%df" % tmp_size
        unpacked = struct.unpack(unpack_format1, decoded)
        for mz, intensity in zip(unpacked[::2], unpacked[1::2]):
            if intensity > self.abs_peak_cutoff:
                self.db_session.add(
                    Peak(scanid=scan.scanid, mz=mz, intensity=intensity))

    def store_mgf(self, mgf_file):
        mgf = open(mgf_file, 'r')
        while mgf.readline()[:10] != 'BEGIN IONS':
            pass
        peaklist = []
        rt=0.0
        while True:
            line = mgf.readline()
            if line[:8] == "PEPMASS=":
                precursormz = float(line.split()[0][8:])
                precursorintensity = int(float(line.split()[1]))
            elif line[:12] == "RTINSECONDS=":
                rt = float(line[12:])/60
            elif line[:8] == 'END IONS':
                break
            else:
                try:
                    mz = float(line.split()[0])
                    intensity = int(float(line.split()[1]))
                    peaklist.append([mz,intensity])
                except:
                    pass
        self.store_peak_list(1, rt, precursormz, [precursormz,precursorintensity], peaklist)

    def store_peak_list(self, scanid, rt, precursormz, basepeak, peaklist):
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
        self.db_session.add(Peak(scanid=scanid, mz=precursormz, intensity=basepeak[1]))
        self.db_session.add(Scan(
            scanid=scanid + 1,
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
            self.db_session.add(Peak(scanid=scanid + 1, mz=peak[0], intensity=peak[1]))
        self.db_session.commit()

    def store_manual_tree(self, manual_tree, tree_type):
        """ Store scans and peaks from mass tree formatted data in scans and peaks tables """
        # tree_type: 0 for mass tree, -1 and 1 for formula trees with negative or positive ionisation mode respectively
        tree_string = open(manual_tree).read()
        tree_string = tree_string.replace(' ', '').\
                                replace('\t', '').\
                                replace('\r', '').\
                                replace('(\n', '(').\
                                replace(',\n', ',').\
                                replace('\n)', ')').\
                                replace('\n', ',')
        tree_list = re.split('([\,\(\)])',tree_string)
        self.global_scanid = 1
        self.store_manual_subtree(tree_list, 0, 0, 0, 1, tree_type)

    def store_manual_subtree(self, tree_list, precursor_scanid, precursor_mz, precursor_intensity, mslevel, tree_type):
        lowmz = None
        highmz = None
        basepeakmz = None
        basepeakintensity = None
        scanid = self.global_scanid
        npeaks = 0
        while len(tree_list) > 0 and tree_list[0] != ')':
            tree_item = tree_list.pop(0)
            if tree_item.find(':') >= 0:
                try:
                    mz, intensity = tree_item.split(':')
                except:
                    raise FileFormatError('Corrupt Tree format ...')
                if tree_type != 0:
                    mz = self.mass_from_formula(mz) - tree_type * pars.elmass
                self.db_session.add(Peak(scanid=scanid, mz=mz, intensity=intensity))
                npeaks += 1
                if lowmz is None or mz < lowmz:
                    lowmz = mz
                if highmz is None or mz > highmz:
                    highmz = mz
                if basepeakintensity is None or intensity > basepeakintensity:
                    basepeakmz = mz
                    basepeakintensity = intensity
            elif tree_item == '(':
                self.global_scanid += 1
                self.store_manual_subtree(tree_list, scanid, mz, intensity, mslevel + 1, tree_type)
            elif tree_item != ',' and tree_item != '':
                raise FileFormatError('Corrupt Tree format ...')
        if npeaks > 0:
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
        if len(tree_list) > 0:
            tree_list.pop(0)

    def mass_from_formula(self, form):
        mass = 0.0
        while len(form) > 0:
            if form[:2] in pars.mims:
                m = pars.mims[form[:2]]
                form = form[2:]
            elif form[:1] in pars.mims:
                m = pars.mims[form[:1]]
                form = form[1:]
            else:
                raise FileFormatError('Element not allowed in formula tree: ' + form)
            x = 0
            while len(form) > x and form[x] in '0123456789':
                x += 1
            if x > 0:
                n = int(form[:x])
            else:
                n = 1
            mass += m * n
            form = form[x:]
        return mass


class AnnotateEngine(object):

    """ Engine to perform MAGMa annotation of the MS/MS data based candidate molecules present
        or retrieved from a chemical database """

    def __init__(self, db_session, skip_fragmentation, max_broken_bonds, max_water_losses,
                 ms_intensity_cutoff, msms_intensity_cutoff, use_all_peaks, adducts=None, max_charge=1, call_back_url=None):
        self.db_session = db_session
        try:
            rundata = self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.skip_fragmentation is None:
            rundata.skip_fragmentation = skip_fragmentation
        if rundata.max_broken_bonds is None:
            rundata.max_broken_bonds = max_broken_bonds
        if rundata.max_water_losses is None:
            rundata.max_water_losses = max_water_losses
        if rundata.ms_intensity_cutoff is None:
            rundata.ms_intensity_cutoff = ms_intensity_cutoff
        if rundata.msms_intensity_cutoff is None:
            rundata.msms_intensity_cutoff = msms_intensity_cutoff
        if rundata.use_all_peaks is None:
            rundata.use_all_peaks = use_all_peaks
        self.db_session.add(rundata)
        self.db_session.commit()
        self.ionisation_mode = rundata.ionisation_mode
        self.skip_fragmentation = rundata.skip_fragmentation
        self.max_broken_bonds = rundata.max_broken_bonds
        self.max_water_losses = rundata.max_water_losses
        self.ms_intensity_cutoff = rundata.ms_intensity_cutoff
        self.msms_intensity_cutoff = rundata.msms_intensity_cutoff
        if rundata.mz_precision is None:
            raise DataProcessingError('No MS data parameters read.')
        self.mz_precision = rundata.mz_precision
        self.precision = 1 + rundata.mz_precision / 1e6
        self.mz_precision_abs = rundata.mz_precision_abs
        self.precursor_mz_precision = rundata.precursor_mz_precision
        self.use_all_peaks = rundata.use_all_peaks

        self.scans = []

        if call_back_url is not None:
            self.call_back_engine = CallBackEngine(call_back_url)
        else:
            self.call_back_engine = None
        logger.info('SETTINGS FOR MATCHING PRECURSOR IONS AND CANDIDATES:')
        logger.info('Ionisation mode: ' + str(self.ionisation_mode))
        if self.ionisation_mode == 1:
            iontypes = ['+H']
        if self.ionisation_mode == -1:
            iontypes = ['-H']
        if adducts is not None:
            for i in adducts.split(','):
                iontypes.append('+' + i)
        self.ions = self.generate_ions(iontypes, max_charge)
        logger.info('Precursor intensity threshold: ' + str(self.ms_intensity_cutoff))
        logger.info('Maximum relative m/z error (ppm): ' + str(self.mz_precision))
        logger.info('Maximum absolute m/z error (Da): ' + str(self.mz_precision_abs) + '\n')

    def generate_ions(self, iontypes, maxcharge):
        """ Generate ions according to allowed number of charges and adducts
            Returns list (one entry per absolute charge) of dictionaries with possible ions {delta mass:symbol}
            Example (negative mode with -a Cl and -m 2):
            [
             {0: '[M]-'},                                           # delta mass 0, (=> candidate molecules already charged)
             {-1.0078250321: '[M-H]-', 34.96885271: '[M+Cl]-'},     # singly charged
             {-2.0156500642: '[M-H-H]2-', 33.9610276779: '[M+Cl-H]2-', 69.93770542: '[M+Cl+Cl]2-'}  #doubly charged, etc.
            ]
            """
        ions = [{0: ''}]
        for c in range(0, maxcharge):
            ions.append({})
            for ionmass in ions[c]:
                for i in iontypes:
                    if i not in pars.ionmasses[self.ionisation_mode]:
                        raise DataProcessingError(
                            'Invalid adduct: ' + i + ' for ionisation mode: ' + str(self.ionisation_mode))
                    ions[c + 1][ionmass + pars.ionmasses[self.ionisation_mode][i]] = ions[c][ionmass] + i
        for c in range(maxcharge + 1):
            for ionmass in ions[c]:
                ions[c][ionmass] = '[M' + ions[c][ionmass] + ']' + str(c) * (c > 1) + \
                    '-' * (self.ionisation_mode < 0) + '+' * (self.ionisation_mode > 0)
        for c in range(1, maxcharge + 1):
            logger.info('Adducts with charge ' + str(self.ionisation_mode * c) + ': ' + str(ions[c].values()))
        return ions

    def build_spectrum(self, dbscan):
        """ Create in memory ScanType object from db Scan object, with list of peaks (of PeakType)
            Return ScanType object and a list of depths of the derived spectral trees """
        scan = types.ScanType(dbscan.scanid, dbscan.mslevel)
        logger.debug('Building scan ' + str(dbscan.scanid))
        if scan.mslevel == 1:
            cutoff = self.ms_intensity_cutoff
        else:
            cutoff = dbscan.basepeakintensity * self.msms_intensity_cutoff / 100
        dbpeaks = self.db_session.query(Peak).filter(
            Peak.scanid == scan.scanid).filter(Peak.intensity >= cutoff).all()
        for dbpeak in dbpeaks:
            scan.peaks.append(types.PeakType(
                dbpeak.mz, dbpeak.intensity, scan.scanid, pars.missingfragmentpenalty * (dbpeak.intensity**0.5)))
        dbchildscans = self.db_session.query(Scan).filter(Scan.precursorscanid == scan.scanid).all()
        max_depth = []  # list of tree depths of all children
        for dbchildscan in dbchildscans:
            # find the highest peak that qualifies as precursor peak for the child spectrum
            prec_intensity = 0.0
            for peak in scan.peaks:
                if peak.intensity > prec_intensity and \
                            -self.precursor_mz_precision < dbchildscan.precursormz - peak.mz < self.precursor_mz_precision:
                    prec_peak = peak
                    prec_intensity = peak.intensity
            # or add the precursor peak as a new peak in the current spectrum
            if prec_intensity == 0.0:
                if dbchildscan.precursorintensity >= cutoff:
                    scan.peaks.append(types.PeakType(dbchildscan.precursormz, dbchildscan.precursorintensity,
                                                     scan.scanid, pars.missingfragmentpenalty * (dbchildscan.precursorintensity**0.5)))
                    prec_peak = scan.peaks[-1]
                else:
                    continue
            max_depth.append(1)
            # process the child spectrum
            prec_peak.childscan, child_depths = self.build_spectrum(dbchildscan)
            for childpeak in prec_peak.childscan.peaks:
                prec_peak.missing_fragment_score += childpeak.missing_fragment_score
            if len(child_depths) > 0:
                max_depth[-1] += max(child_depths)
        return scan, max_depth

    def build_spectra(self, scans='all'):
        """ Build list of scans (of ScanType) requested for annotation
            Precursor peaks are indexed based on integer mass """
        logger.info('BUILDING SPECTRAL TREES')
        ndepths = {}
        if scans == 'all':
            queryscans = self.db_session.query(Scan).filter(Scan.mslevel == 1).all()
        else:
            queryscans = self.db_session.query(Scan).filter(Scan.mslevel == 1).filter(Scan.scanid.in_(scans)).all()
        for dbscan in queryscans:
            spectrum, max_depth = self.build_spectrum(dbscan)
            for depth in max_depth:
                depth += 1
                if depth > 1:
                    if depth in ndepths:
                        ndepths[depth] += 1
                    else:
                        ndepths[depth] = 1
            self.scans.append(spectrum)
        logger.info(str(len(self.scans)) + ' MS1 spectra')
        for depth in ndepths:
            logger.info(str(ndepths[depth]) + ' spectral trees of depth ' + str(depth))
        logger.info('')
        self.indexed_peaks = {}   # sets of peaks for each integer m/z value
        for scan in self.scans:
            for peak in scan.peaks:
                if not ((not self.use_all_peaks) and peak.childscan is None):
                    int_mass = int(round(peak.mz))
                    if int_mass not in self.indexed_peaks:
                        self.indexed_peaks[int_mass] = set([])
                    self.indexed_peaks[int_mass].add(peak)
                    # write mass_tree for selected level 1 peak
                    logger.debug('Mass_tree for scan ' + str(scan.scanid) + ', m/z=' + str(peak.mz) + ':\n' +
                                 self.write_mass_tree(peak))

    def write_mass_tree(self, peak):
        peak_string = "%.6f: %i" % (peak.mz, peak.intensity)
        if peak.childscan is not None:
            peak_string += ' ('
            n = 0
            for childpeak in peak.childscan.peaks:
                if n > 0:
                    peak_string += ", "
                peak_string += self.write_mass_tree(childpeak)
                n += 1
            peak_string += ')'
        return peak_string

    def get_db_candidates(self, query_engine, max_mim=""):
        """ Query given query engine to retrieve all relevant candidate molecules,
            store in molecules table and return corresponding molids """
        logger.info('RETRIEVING CANDIDATE MOLECULES FROM: ' + str(query_engine.name))
        if max_mim == '':
            max_mim = '1200'
        logger.info('Mass limit: ' + str(max_mim))

        struct_engine = StructureEngine(self.db_session)

        # build sorted list of query masses
        mzs = []
        for scan in self.scans:
            for peak in scan.peaks:
                if not ((not self.use_all_peaks) and peak.childscan is None):
                    mzs.append(peak.mz)
        mzs.sort()
        # build non-overlapping set of queries around these masses
        candidate = {}
        for mz in mzs:
            for charge in range(1, len(self.ions)):
                for ionmass in self.ions[charge]:
                    ql = int(1e6 * ((min(mz / self.precision, mz - self.mz_precision_abs) - ionmass) *
                                    charge + self.ionisation_mode * pars.elmass))
                    qh = int(1e6 * ((max(mz * self.precision, mz + self.mz_precision_abs) - ionmass) *
                                    charge + self.ionisation_mode * pars.elmass))
                    result = query_engine.query_on_mim(ql, qh, 0)
                    for molecule in result:
                        # remove duplicates
                        candidate[molecule.inchikey14] = molecule
                    logger.debug(str(ql) + ',' + str(qh) + ' --> ' + str(len(candidate)) + ' candidates')
                # include singly charged candidates from database
                for ionmass in self.ions[charge - 1]:
                    ql = int(1e6 * ((min(mz / self.precision, mz - self.mz_precision_abs) - ionmass) *
                                    charge + self.ionisation_mode * pars.elmass))
                    qh = int(1e6 * ((max(mz * self.precision, mz + self.mz_precision_abs) - ionmass) *
                                    charge + self.ionisation_mode * pars.elmass))
                    result = query_engine.query_on_mim(ql, qh, self.ionisation_mode)
                    for molecule in result:
                        # remove duplicates
                        candidate[molecule.inchikey14] = molecule
                    logger.debug(str(ql) + ',' + str(qh) + ' --> ' + str(len(candidate)) + ' candidates')
        logger.info('Storing ' + str(len(candidate)) + ' candidates')
        # in case of an empty database, no check for existing duplicates needed
        check_duplicates = (self.db_session.query(Molecule.molid).count() > 0)
        logger.debug('check duplicates: ' + str(check_duplicates))
        molids = set([])
        for molecule in candidate.itervalues():
            molid = struct_engine.add_molecule(molecule, check_duplicates=check_duplicates, merge=True)
            molids.add(molid)

        # All candidates are stored in dbsession, resulting molids are returned
        self.db_session.commit()
        return molids

    def search_structures(self, molids=None, ncpus=1, fast=False, time_limit=None):
        """ Match candidate molecules with precursor ions, find substructures
            for fragment peaks and calculate candidate scores """
        logger.info('MATCHING CANDIDATE MOLECULES')
        global fragid
        fragid = self.db_session.query(func.max(Fragment.fragid)).scalar()
        if fragid is None:
            fragid = 0
        ppservers = ()
        logger.info('calculating on ' + str(ncpus) + ' cpus')
        job_server = pp.Server(ncpus, ppservers=ppservers)
        if molids is None:
            if time_limit is None:
                metabdata = self.db_session.query(Molecule.molid).order_by(desc(Molecule.molid)).all()
            else:
                metabdata = self.db_session.query(Molecule.molid).order_by(Molecule.refscore).all()
            molids = [x[0] for x in metabdata]
        # annotate molids in chunks of 500 to avoid errors in db_session.query
        # and memory problems during parallel processing
        total_frags = 0
        total_molids = len(molids)
        start_time = time.time()
        while len(molids) > 0:
            ids = set([])
            while len(ids) < 500 and len(molids) > 0:
                ids.add(molids.pop())
            structures = self.db_session.query(Molecule).filter(Molecule.molid.in_(ids)).all()
            jobs = []
            for structure in structures:
                # skip molecule if it has already been used for annotation
                if self.db_session.query(Fragment.fragid).filter(Fragment.molid == structure.molid).count() > 0:
                    logger.warn('Molecule ' + str(structure.molid) + ': Already annotated, skipped')
                    continue
                # collect all peaks with masses within 3 Da range
                molcharge = 0
                # derive charge from molecular formula
                molcharge += 1 * ((structure.formula[-1] == '-' and self.ionisation_mode == -1) or
                                  (structure.formula[-1] == '+' and self.ionisation_mode == 1))
                peaks = set([])
                # select a subset of precursor peaks potentially matching the candidate structure
                for charge in range(1, len(self.ions)):
                    for ionmass in self.ions[charge - molcharge]:
                        int_mass = int(round((structure.mim + ionmass) / charge))
                        try:
                            peaks = peaks.union(self.indexed_peaks[int_mass])
                        except:
                            pass
                        try:
                            peaks = peaks.union(self.indexed_peaks[int_mass - 1])
                        except:
                            pass
                        try:
                            peaks = peaks.union(self.indexed_peaks[int_mass + 1])
                        except:
                            pass
                if len(peaks) == 0:
                    logger.debug('Molecule ' + str(structure.molid) + ': No match')
                    continue
                if fast and structure.natoms <= 64:
                    fragmentation_module = 'magma.fragmentation_cy'
                else:
                    fragmentation_module = 'magma.fragmentation_py'
                # distribute search_structure tasks over different cpus
                jobs.append((structure,
                             job_server.submit(search_structure, (structure.mol,
                                                                  structure.mim,
                                                                  molcharge,
                                                                  peaks,
                                                                  self.max_broken_bonds,
                                                                  self.max_water_losses,
                                                                  self.precision,
                                                                  self.mz_precision_abs,
                                                                  self.use_all_peaks,
                                                                  self.ionisation_mode,
                                                                  self.skip_fragmentation,
                                                                  (fast and structure.natoms <=
                                                                   64),
                                                                  self.ions
                                                                  ), (), (
                                 "magma.types",
                                 "magma.pars",
                                 fragmentation_module
                             )
                             )))
            # process results from search_structure tasks
            count = 0
            for structure, job in jobs:
                raw_result = job(raw_result=True)
                result, sout = pickle.loads(raw_result)
                logger.debug(sout)
                (hits, frags) = result
                total_frags += frags
                structure.nhits = len(hits)
                self.db_session.add(structure)
                if len(hits) == 0:
                    logger.debug('Molecule ' + str(structure.molid) + ': No match')
                else:
                    logger.debug('Molecule ' + str(structure.molid) + ': ' +
                                 structure.name.encode('utf-8') + ' -> ' + str(frags) + ' fragments')
                    for hit in hits:
                        score = self.store_hit(hit, structure.molid, 0)
                        logger.debug('Scan: ' + str(hit.scan) + ' - Mz: ' + str(hit.mz) + ' - ' + 'Score: ' + str(score))
                self.db_session.flush()
                count += 1
                elapsed_time = time.time() - start_time
                if self.call_back_engine is not None:
                    status = 'Annotation: %d / %d candidate molecules processed  (%d%%)' % (total_molids - len(molids) - len(ids) + count,
                        total_molids, 100.0 * (total_molids - len(molids) - len(ids) + count) / total_molids)
                    self.call_back_engine.update_callback_url(status, elapsed_time, time_limit)
                if time_limit and elapsed_time > time_limit * 60:
                    # break out of while-loop
                    molids = []
                    if self.call_back_engine is not None:
                        self.call_back_engine.update_callback_url(
                            'Annotation stopped: time limit exceeded', force=True)
                    logger.warn('Annotation stopped: time limit exceeded')
                    break
            self.db_session.commit()
            logger.info(str(total_molids - len(molids)) + ' molecules processed')
        if self.call_back_engine is not None:
            self.call_back_engine.update_callback_url(
                'Annotation completed', force=True)
        logger.info(str(total_frags) + ' fragments generated in total.')
        nmols = (self.db_session.query(Fragment.molid).filter(Fragment.parentfragid == 0).distinct().count())
        nprecursors = (self.db_session.query(Fragment.scanid, Fragment.mz).filter(Fragment.parentfragid == 0).distinct().count())
        logger.info(str(nmols) + ' Molecules matched with ' + str(nprecursors) + ' precursor ions, in total\n')
        job_server.destroy()

    def store_hit(self, hit, molid, parentfragid):
        """ Store candidate molecule and its substructures in fragments table """
        global fragid
        fragid += 1
        currentFragid = fragid
        score = hit.score
        deltappm = None
        if score is not None:
            score = score / hit.intensity_weight
            charge = 1
            if hit.ion[-2] in '123456789':
                charge = int(hit.ion[-2])  # TODO store charge of a hit explicitly
            deltappm = (hit.mz - (hit.mass + hit.deltaH) / charge +
                        self.ionisation_mode * pars.elmass) / hit.mz * 1e6
        self.db_session.add(Fragment(
            molid=molid,
            scanid=hit.scan,
            mz=hit.mz,
            mass=hit.mass,
            score=score,
            parentfragid=parentfragid,
            atoms=unicode(hit.atomstring),
            smiles=unicode(hit.smiles),
            deltah=hit.deltaH,
            deltappm=deltappm,
            formula=unicode(hit.formula+'<br>'+hit.ion)
            ))
        if len(hit.besthits) > 0:
            for childhit in hit.besthits:
                if childhit is not None:
                    self.store_hit(childhit, molid, currentFragid)
        return score


class PubChemEngine(object):

    """Engine to retrieve candidate molecules from PubChem"""

    def __init__(self, db, dbfilename='', max_64atoms=False, incl_halo='', min_refscore='', online=True):
        self.name = db
        self.incl_halo = False
        if incl_halo != '' and incl_halo != False:
            self.incl_halo = True
        if config.getboolean('magma job', 'structure_database.online') and online:
            self.query = self.query_online
            self.service = config.get('magma job', 'structure_database.service')+'/'+db
        else:
            self.query = self.query_local
            databases = {'pubchem': 'structure_database.pubchem', 'kegg': 'structure_database.kegg'}
            halo_databases = {'pubchem': 'structure_database.pubchem_halo', 'kegg': 'structure_database.kegg_halo'}
            if dbfilename == '':
                dbfilename = config.get('magma job', databases[db])
            self.conn = sqlite3.connect(dbfilename)
            self.conn.text_factory = str
            self.c = self.conn.cursor()
            if self.incl_halo == True:
                halo_filename = config.get('magma job', halo_databases[db])
                self.connh = sqlite3.connect(halo_filename)
                self.connh.text_factory = str
                self.ch = self.connh.cursor()
            self.where = ''
            if min_refscore != '':
                self.where += ' AND refscore >= ' + min_refscore
            if max_64atoms:
                self.where += ' AND natoms <= 64'

    def query_online(self, low, high, charge):
        r = requests.post(self.service, data=json.dumps([low, high, charge, self.incl_halo]))
        try:
            result = r.json()
        except:
            result = r.json # for compatibility with older version of requests
        return result

    def query_local(self, low, high, charge):
        """ Return all molecules with given charge from HMDB between low and high mass limits """
        result = self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where,
                                (charge, low, high)).fetchall()
        if self.incl_halo:
            result += self.ch.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where,
                                (charge, low, high)).fetchall()
        return result

    def query_on_mim(self, low, high, charge):
        """ Return all molecules with given charge between low and high mass limits """
        molecules = []
        result = self.query(low, high, charge)
        for (cid, mim, charge, natoms, molblock, inchikey, smiles, molform, name, refs, logp) in result:
            if self.name == 'pubchem':
                refscore = refs
                reference = '<a href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&cmd=Link&'+\
                            'LinkName=pccompound_pccompound_sameisotopic_pulldown&from_uid='+str(cid)+'" target="_blank">'+str(cid)+' (PubChem)</a>'
            elif self.name == 'kegg':
                refscore = None
                keggids = refs.split(',')
                reference = '<a href="http://www.genome.jp/dbget-bin/www_bget?cpd:' + \
                    keggids[0] + '" target="_blank">' + keggids[0] + ' (Kegg)</a>'
                for keggid in keggids[1:]:
                    reference += '<br><a href="http://www.genome.jp/dbget-bin/www_bget?cpd:' + \
                        keggid + '" target="_blank">' + keggid + ' (Kegg)</a>'
            else:
                exit(self.name)
            molecules.append(get_molecule(
                           zlib.decompress(base64.decodestring(molblock)),
                           name+' (' + str(cid) + ')',
                           refscore,
                           0,
                           mim=float(mim / 1e6),
                           natoms=natoms,
                           molform=molform,
                           inchikey14=inchikey,
                           smiles=smiles,
                           reference=reference,
                           logp=float(logp) / 10.0,
                           ))
        return molecules

    def check_inchi(self, mim, inchikey14):
        """ Function to look up uploaded structures in PubChem based on inchikey. Returns refscore and reference.
            Only available if PubChem database is installed locally""" 
        if config.getboolean('magma job', 'structure_database.online'):
            return False
        self.c.execute('SELECT cid,name,refscore FROM molecules WHERE charge IN (-1,0,1) AND mim between ? and ? and inchikey = ?',
                       (int(mim * 1e6) - 1, int(mim * 1e6) + 1, inchikey14))
        result = self.c.fetchall()
        if len(result) > 0:
            cid, name, refscore = result[0]
            reference = '<a href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&cmd=Link&LinkName=pccompound_pccompound_sameisotopic_pulldown&from_uid=' +\
                str(cid) + '" target="_blank">' + str(cid) + ' (PubChem)</a>'
            return [name, reference, refscore]
        else:
            return False


class HmdbEngine(object):

    """Engine to retrieve candidate molecules from HMDB"""

    def __init__(self, dbfilename='', max_64atoms=False, online=True):
        self.name = 'Human Metabolite Database'
        if config.getboolean('magma job', 'structure_database.online') and online:
            self.query = self.query_online
            self.service = config.get('magma job', 'structure_database.service')+'/hmdb'
        else:
            self.query = self.query_local
            if dbfilename == '':
                dbfilename = config.get('magma job', 'structure_database.hmdb')
            self.where = ''
            if max_64atoms == True:
                self.where += ' AND natoms <= 64'
            self.conn = sqlite3.connect(dbfilename)
            self.conn.text_factory = str
            self.c = self.conn.cursor()
        
    def query_online(self, low, high, charge):
        r = requests.post(self.service, data=json.dumps([low, high, charge, True]))
        try:
            result = r.json()
        except:
            result = r.json # for compatibility with older version of requests
        return result

    def query_local(self, low, high, charge):
        """ Return all molecules with given charge from HMDB between low and high mass limits """
        return self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where,
                                (charge, low, high)).fetchall()

    def query_on_mim(self, low, high, charge):
        molecules = []
        result = self.query(low, high, charge)
        for (cid, mim, charge, natoms, molblock, inchikey, smiles, molform, name, reference, logp) in result:
            hmdb_ids = reference.split(',')
            hmdb_refs = '<a href="http://www.hmdb.ca/metabolites/' + \
                hmdb_ids[0] + '" target="_blank">' + hmdb_ids[0] + ' (HMDB)</a>'
            for hmdb_id in hmdb_ids[1:]:
                hmdb_refs += '<br><a href="http://www.hmdb.ca/metabolites/' + \
                    hmdb_id + '" target="_blank">' + hmdb_id + ' (HMDB)</a>'
            molecules.append(get_molecule(
                           zlib.decompress(base64.decodestring(molblock)),
                           name + ' (' + str(cid) + ')',
                           None,
                           0,
                           mim=float(mim / 1e6),
                           natoms=natoms,
                           molform=molform,
                           inchikey14=inchikey,
                           smiles=smiles,
                           reference=hmdb_refs,
                           logp=float(logp) / 10.0,
                           ))
        return molecules


class ExportMoleculesEngine(object):

    def __init__(self, db_session):
        self.db_session = db_session

    def export_molecules(self, output_format='sdf', filename=None, columns=None, sortcolumn='refscore', descend=True):
        """ Write SDFile or smiles with candidate molecules to filename (or stdout).
            If data is for a single percursor ion: also provide candidate scores and sort accordingly """
        if filename is None:
            file = sys.stdout
        else:
            file = open(filename, 'w')
        # If data is for a single percursor ion: also provide candidate scores and sort accordingly
        nprecursors = self.db_session.query(Fragment.mz, Fragment.scanid).\
                      filter(Fragment.parentfragid == 0).distinct().count()
        if nprecursors == 1:
            result = self.db_session.query(Molecule, Fragment.score).\
                     filter(Molecule.molid == Fragment.molid).\
                     filter(Fragment.parentfragid == 0).\
                     order_by(Fragment.score, desc(Molecule.refscore)).all()
        else:
            if descend:
                result = self.db_session.query(Molecule, Molecule.molid).order_by(desc(sortcolumn)).all()
            else:
                result = self.db_session.query(Molecule, Molecule.molid).order_by(sortcolumn).all()
        for molecule, value in result:
            if output_format == 'sdf':
                file.write(molecule.mol)
                if nprecursors == 1:
                    file.write('> <score>\n%.5f\n\n' % value)
                if columns is None:
                    columns = dir(molecule)
                for column in columns:
                    if column[:1] != '_' and column != 'mol' and column != 'metadata' and column != 'fragments':
                        file.write('> <' + column + '>\n' + str(molecule.__getattribute__(column)) + '\n\n')
                file.write('$$$$\n')
            else:
                file.write(molecule.smiles)
                if nprecursors == 1:
                    file.write(' score=%.5f' % value)
                if columns is None:
                    columns = ['name','refscore','formula','mim']
                for column in columns:
                    if column[:1] != '_' and column != 'mol' and column != 'metadata' and column != 'fragments' and column != 'smiles':
                        file.write(' ' + column + '=' + str(molecule.__getattribute__(column)).replace(" ","_"))
                file.write('\n')
                
        file.close()

    def export_assigned_molecules(self, filename=None):
        """ Write SDFile with assigned candidate molecules to filename (or stdout).
            Include scanid, RT and m/z of assigned peaks """
        if filename is None:
            file = sys.stdout
        else:
            file = open(filename, 'w')
        result = self.db_session.query(Molecule, Peak, Scan.rt).\
                 filter(Molecule.molid == Peak.assigned_molid).filter(Peak.scanid == Scan.scanid).all()
        for molecule, peak, rt in result:
            file.write(molecule.mol)
            cm = dir(molecule)
            cp = dir(peak)
            for column in cm:
                if column[:1] != '_' and column != 'mol' and column != 'metadata' and column != 'fragments':
                    file.write('> <' + column + '>\n' + str(molecule.__getattribute__(column)) + '\n\n')
            for column in cp:
                if column[:1] != '_' and column != 'metadata' and column != 'scan':
                    file.write('> <' + column + '>\n' + str(peak.__getattribute__(column)) + '\n\n')
            file.write('> <rt>\n' + str(rt) + '\n\n')
            file.write('$$$$\n')
        file.close()


def search_structure(mol, mim, molcharge, peaks, max_broken_bonds, max_water_losses, precision,
                     mz_precision_abs, use_all_peaks, ionisation_mode, skip_fragmentation, fast, ions):
    """ Match a candidate molecule with precursor ions.
        Return a list of hits (=hierarchical trees of (sub)structures and scores) """
    pars = magma.pars
    if fast:
        Fragmentation = magma.fragmentation_cy
    else:
        Fragmentation = magma.fragmentation_py

    def massmatch(peak, mim, molcharge):
        lowmz = min(peak.mz / precision, peak.mz - mz_precision_abs)
        highmz = max(peak.mz * precision, peak.mz + mz_precision_abs)
        # lowmz=peak.mz/precision
        # highmz=peak.mz*precision
        for charge in range(1, len(ions)):
            for ionmass in ions[charge - molcharge]:
                if lowmz <= (mim + ionmass) / charge - ionisation_mode * pars.elmass <= highmz:
                    return [ionmass, ions[charge - molcharge][ionmass]]
        return False

    def gethit(peak, fragment, score, bondbreaks, mass, ionmass, ion):
        try:
            hit = types.HitType(peak, fragment, score, bondbreaks, mass, ionmass, ion)
        except:
            hit = magma.types.HitType(peak, fragment, score, bondbreaks, mass, ionmass, ion)
        # fragment=0 means it is a missing fragment
        if fragment > 0 and peak.childscan is not None and len(peak.childscan.peaks) > 0:
            n_child_peaks = len(peak.childscan.peaks)
            total_score = 0.0
            total_count = 0.0
            for childpeak in peak.childscan.peaks:
                besthit = gethit(childpeak, 0, None, 0, 0, 0, '')
                # calculate m/z value of the neutral form of the fragment (mass of electron added/removed) for matching theoretical masses
                mz_neutral = childpeak.mz + ionisation_mode * pars.elmass
                for childfrag, childscore, childbbreaks, childmass, childH in \
                        fragment_engine.find_fragments(mz_neutral, fragment, precision, mz_precision_abs):
                    if childfrag & fragment == childfrag:
                        ion = '[X' + '+' * (childH > 0) + '-' * (childH < 0) + str(abs(childH)) * (not -2 < childH < 2) + \
                                'H' * (childH != 0) + ']' + '+' * (ionisation_mode > 0) + '-' * (ionisation_mode < 0)
                        childhit = gethit(childpeak, childfrag, childscore * (childpeak.intensity**0.5),
                                          childbbreaks, childmass, childH * pars.Hmass, ion)
                        if besthit.score is None or besthit.score > childhit.score or \
                                (besthit.score == childhit.score and abs(besthit.deltaH) > abs(childhit.deltaH)) or \
                                fragment_engine.score_fragment_rel2parent(besthit.fragment, fragment) > \
                                        fragment_engine.score_fragment_rel2parent(childhit.fragment, fragment):
                            besthit = childhit
                if besthit.score is None:
                    total_score += childpeak.missing_fragment_score
                    # total_score+=missingfragmentpenalty*weight
                else:
                    hit.besthits.append(besthit)
                    total_score += min(besthit.score,
                                       childpeak.missing_fragment_score)
                    # total_score+=min(besthit.score,missingfragmentpenalty)*weight
            hit.score = hit.score + total_score
        return hit

    def add_fragment_data_to_hit(hit):
        if hit.fragment != 0:
            hit.atomstring, hit.atomlist, hit.formula, hit.smiles = fragment_engine.get_fragment_info(
                hit.fragment, hit.deltaH)
            # except:
            #    exit('failed inchi for: '+atomstring+'--'+str(hit.fragment))
            if len(hit.besthits) > 0:
                for childhit in hit.besthits:
                    if childhit is not None:
                        add_fragment_data_to_hit(childhit)

    # main loop
    Fragmented = False
    hits = []
    frags = 0
    for peak in peaks:
        if not ((not use_all_peaks) and peak.childscan is None):
            i = massmatch(peak, mim, molcharge)
            if i != False:
                if not Fragmented:
                    fragment_engine = Fragmentation.FragmentEngine(
                        mol, max_broken_bonds, max_water_losses, ionisation_mode, skip_fragmentation, molcharge)
                    if fragment_engine.accepted():
                        frags = fragment_engine.generate_fragments()
                    Fragmented = True
                if fragment_engine.accepted():
                    hit = gethit(peak, (1 << fragment_engine.get_natoms()) - 1, 0, 0, mim, i[0], i[1])
                    add_fragment_data_to_hit(hit)
                    hits.append(hit)
    return (hits, frags)
