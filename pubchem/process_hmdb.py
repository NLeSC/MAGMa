#!/usr/bin/env python
"Update local HMDB database"

import os
import zlib
import sqlite3
import sys
import argparse
import urllib2
import zipfile
import StringIO
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Descriptors


def version():
    return '1.0'


def process_hmdb(args):
    conn = sqlite3.connect(args.database_dir + '/HMDB_MAGMa.db')
    c = conn.cursor()
    try:
        c.execute("""CREATE TABLE molecules (id TEXT PRIMARY KEY,
                                             mim INTEGER NOT NULL,
                                             charge INTEGER NOT NULL,
                                             natoms INTEGER NOT NULL,
                                             molblock BLOB,
                                             inchikey TEXT,
                                             smiles TEXT,
                                             molform TEXT,
                                             name TEXT,
                                             reference TEXT,
                                             logp INT)""")
        conn.commit()
        print ("HMDB_MAGMa.db created")
    except:
        print ("HMDB_MAGMa.db already exists (or error creating it)")
        exit()

    if args.data_dir == None:
        zf = urllib2.urlopen('http://www.hmdb.ca/downloads/structures.zip')
    else:
        zf = open(args.data_dir + 'structures.zip')
    sdfile = zipfile.ZipFile(StringIO.StringIO(zf.read())).open('structures.sdf')

    memstore = {}
    line = '$$$$'
    while line != "":
        record = []
        amap = {}
        skip = False
        ionized = 0
        # read heading:
        for x in range(4):
            line = sdfile.readline()
            record.append(line)
        if line == "":
            continue
        natoms = int(record[-1][:3])
        nbonds = int(record[-1][3:6])
        bonds = 0
        y = 0
        for x in range(natoms):
            line = sdfile.readline()
            if line[31:33] == 'H ':
                # skip hydrogens
                continue
            y += 1
            amap[x + 1] = y
            if line[31:33] not in ['C ', 'N ', 'O ', 'P ', 'S ', 'F ', 'Cl', 'Br', 'I ']:
                # filter non-organic compounds
                skip = True
            elif line[50:51] != '0':
                # this flag has something to do with polymeric structures
                # and resulted in deviation between calculated and given inchikeys, skip
                skip = True
            elif line[38:39] == '4':
                # radical, resulted in deviation between calculated and given inchikeys
                skip = True
            record.append(line[:42] + '\n')
        for x in range(nbonds):
            line = sdfile.readline()
            a1 = int(line[:3])
            a2 = int(line[3:6])
            # skip bonds involving hydrogens
            if a1 in amap and a2 in amap:
                bonds += 1
                # use bonds with stereoflags set to zero
                record.append('%3i%3i%s  0\n' %
                              (amap[a1], amap[a2], line[6:9]))
        while line != 'M  END\n' and line != '':
            line = sdfile.readline()
            record.append(line)
            if line[:6] == 'M  ISO':
                skip = True
                print 'Skipped isotopically labeled:', record[0][:-1]
        while line != "$$$$\n" and line != "":
            line = sdfile.readline()
            if line == "> <HMDB_ID>\n":
                hmdb_id = str(sdfile.readline()[:-1])
            elif line == "> <GENERIC_NAME>\n":
                molname = str(sdfile.readline()[:-1])
            elif line == "> <INCHI_KEY>\n":
                inchi_key = sdfile.readline()[9:-1]
        if line != "" and skip == False:
            record[3] = repr(y).rjust(3) + repr(bonds).rjust(3) + record[3][6:]
            molblock = ''.join(record)
            mol = Chem.MolFromMolBlock(molblock)
            if mol == None or mol.GetNumAtoms() == 0:
                continue
            smiles = Chem.MolToSmiles(mol)
            if len(Chem.GetMolFrags(mol)) > 1:
                print 'complex:', hmdb_id, smiles
                continue
            conf = mol.GetConformer(0)
            molblock = zlib.compress(''.join(record))
            molform = Chem.rdMolDescriptors.CalcMolFormula(mol)
            mim = Chem.rdMolDescriptors.CalcExactMolWt(mol)
            charge = 0
            if '-' in molform:
                if molform[-1] == '-':
                    charge = -1
                else:
                    continue
            elif '+' in molform:
                if molform[-1] == '+':
                    charge = 1
                else:
                    continue
            if mim > 1200.0:
                print 'molecule to heavy:', hmdb_id, smiles
                continue
            natoms = mol.GetNumHeavyAtoms()
            logp = Chem.Crippen.MolLogP(mol)
            inchikey = Chem.AllChem.InchiToInchiKey(
                AllChem.MolToInchi(mol))[:14]
            if inchikey != inchi_key[:14]:
                print 'given inchikey does not match calculated inchikey, skipped:', hmdb_id, smiles
                continue
            ionized = 0
            for x in ['C(=O)[O-]', '[NH+]', '[NH2+]', '[NH3+]', '[NH4+]']:
                if smiles.find(x) >= 0:
                    ionized = 1
            if inchikey in memstore:
                dbid, reference, dbionized = memstore[inchikey]
                reference = reference + ',' + hmdb_id
                print 'Duplicates:', reference, molname
                if dbionized > ionized:  # prefer non-ionized CID's
                    c.execute('''UPDATE molecules SET id=?, mim=?, charge=?, molblock=?, smiles=?,
                                 molform=?, name=?, reference=?, logp=? WHERE id == ?''', (
                                    hmdb_id,
                                    int(mim * 1e6),
                                    charge,
                                    buffer(molblock),
                                    unicode(smiles),
                                    unicode(molform),
                                    unicode(molname, 'utf-8', 'xmlcharrefreplace'),
                                    unicode(reference),
                                    int(logp * 10),
                                    dbid))
                    memstore[inchikey] = (hmdb_id, reference, ionized)
                else:
                    c.execute('UPDATE molecules SET reference=? WHERE id == ?', (
                        unicode(reference),dbid))
                    memstore[inchikey] = (dbid, reference, dbionized)
            else:
                c.execute('''INSERT INTO molecules (id, mim, charge, natoms, molblock, inchikey,
                             smiles,molform,name,reference,logp) VALUES (?,?,?,?,?,?,?,?,?,?,?)''', (
                                hmdb_id,
                                int(mim * 1e6),
                                charge,
                                int(natoms),
                                buffer(molblock),
                                unicode(inchikey),
                                unicode(smiles),
                                unicode(molform),
                                unicode(molname, 'utf-8', 'xmlcharrefreplace'),
                                unicode(hmdb_id),
                                int(logp * 10)))
                memstore[inchikey] = (hmdb_id, hmdb_id, ionized)
    conn.commit()

    print "Creating index ..."
    c.execute('PRAGMA temp_store = 2')
    c.execute(
        'CREATE INDEX idx_cover ON molecules (charge,mim,natoms,reference,molform,inchikey,smiles,name,molblock,logp)')
    conn.commit()

# main
mainparser = argparse.ArgumentParser(description=__doc__)
mainparser.add_argument(
    '-v', '--version', action='version', version='%(prog)s ' + version())
mainparser.add_argument(
    '-d', '--data_dir', help="""Directory where HMDB structures.zip file is stored
                                (default: %(default)s, attempt to read directly from HMDB server)""", default=None, type=str)
mainparser.add_argument(
    '-b', '--database_dir', help="Directory where HMDB database is stored (default: %(default)s)", default="./", type=str)
mainparser.set_defaults(func=process_hmdb)

args = mainparser.parse_args(sys.argv[1:])
args.func(args)
