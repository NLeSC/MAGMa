#!/usr/bin/env python
"Generate local PubChem compound databases for MAGMa"

import os
import re
import zlib
import gzip
import sqlite3
import time
import sys
import argparse
import base64


def version():
    return '1.0'


def process(args):
    """Update local PubChem database.
    This has to be run twice, with and without -f option,
    to generate db's for halogenated and non-halogenated compounds respectively"""
    kegg_info = {}
    if args.kegg != None:
        kegg_info = parse_kegg_file(args.kegg)
        print len(kegg_info), 'kegg compound names read'
    if not args.skip_names:
        create_names_db(args.data_dir, args.database_dir)
    create_pubchem_dbs(args.data_dir, args.database_dir, kegg_info, args.halogens)


def create_names_db(data_dir, dbs_dir):
    "Generate Pubchem_names.db which is used later by create_pubchem_dbs"
    try:
        os.remove(dbs_dir + "/Pubchem_Names.db")
    except:
        pass
    conn = sqlite3.connect(dbs_dir + "/Pubchem_Names.db")
    curs = conn.cursor()
    curs.execute("CREATE TABLE names (cid INTEGER PRIMARY KEY, name TEXT, refscore INTEGER)")
    meshfile = gzip.open(data_dir + "/ftp.ebi.ac.uk/pub/databases/pubchem/Compound/Extras/CID-MeSH.gz")
    mesh = {}
    for line in meshfile:
        splitline = line.split('\t')
        mesh[splitline[0]] = splitline[1][:-1]
    meshfile.close()
    namesfile = gzip.open(
        data_dir + '/ftp.ebi.ac.uk/pub/databases/pubchem/Compound/Extras/CID-Synonym-filtered.gz')
    sidfile = gzip.open(
        data_dir + '/ftp.ebi.ac.uk/pub/databases/pubchem/Compound/Extras/CID-SID.gz')
    splitsid = sidfile.readline().split("\t")
    splitnames = namesfile.readline().split("\t")
    while splitsid != ['']:
        curr_cid = splitsid[0]
        sid_count = 1
        name = ''
        splitsid = sidfile.readline().split("\t")
        while splitsid[0] == curr_cid:
            sid_count += 1
            splitsid = sidfile.readline().split("\t")
        if splitnames[0] == curr_cid:
            name = splitnames[1][:-1]
        if curr_cid in mesh:
            name = mesh[curr_cid]
        while splitnames[0] == curr_cid:
            splitnames = namesfile.readline().split("\t")
        try:
            curs.execute('INSERT INTO names (cid, name, refscore) VALUES (?,?,?)', (int(
                curr_cid), unicode(name, 'utf-8', 'xmlcharrefreplace'), sid_count))
        except:
            print 'Compound name:\n', name, '\ncould not be encoded.'
    conn.commit()


def update(args):
    """Update downloads from PubChem server"""
    command1 = ""
    if args.data_dir != None:
        command1 = "cd " + args.data_dir + '; '
    command2 = 'wget --mirror --accept "*.sdf.gz" ftp://ftp.ebi.ac.uk/pub/databases/pubchem/Compound/Extras/CID-MeSH.gz'
    os.system(command1 + command2)
    command2 = 'wget --mirror --accept "*.sdf.gz" ftp://ftp.ebi.ac.uk/pub/databases/pubchem/Compound/Extras/CID-Synonym-filtered.gz'
    os.system(command1 + command2)
    command2 = 'wget --mirror --accept "*.sdf.gz" ftp://ftp.ebi.ac.uk/pub/databases/pubchem/Compound/Extras/CID-SID.gz'
    os.system(command1 + command2)
    command2 = 'wget --mirror --accept "*.sdf.gz" ftp://ftp.ebi.ac.uk/pub/databases/pubchem/Compound/CURRENT-Full/SDF/'
    os.system(command1 + command2)


def parse_kegg_file(kegg_file):
    """Generate kegg_info dictionary from kegg_file
       kegg_info[cid]=(keggname,keggid)"""
    kf = open(kegg_file)
    kegg_info = {}
    line = kf.readline()
    while line != "":
        line = kf.readline()
        keggname = line.split(';')[0].split('. ')[1]
        line = kf.readline()
        line = kf.readline()
        if line.find('CID') >= 0:
            cid = int(line.split()[3])
            keggid = line.split()[5]
            if keggid[0] == 'C':
                kegg_info[cid] = (keggname, keggid)
        line = kf.readline()
    return kegg_info


def create_pubchem_dbs(data_dir, dbs_dir, kegg_info, halogens):
    """Generate Pubchem_MAGMa.db and Pubchem_MAGMa_kegg.db in dbs_dir from data in data_dir"""

    pubchem_dir = data_dir + \
        "/ftp.ebi.ac.uk/pub/databases/pubchem/Compound/CURRENT-Full/SDF/"
    # open database with compound names
    conn_names = sqlite3.connect(dbs_dir + "/Pubchem_Names.db")
    curs_names = conn_names.cursor()

    # generate a database with a list of all compound SDFs
    try:
        os.remove(dbs_dir + "/Pubchem_Listing.db")
    except:
        pass
    conn_list = sqlite3.connect(dbs_dir + '/Pubchem_Listing.db')
    curs_list = conn_list.cursor()
    curs_list.execute(
        'CREATE TABLE files (file TEXT UNIQUE NOT NULL, processed INTEGER, time REAL)')
    conn_list.commit()
    listing = os.listdir(pubchem_dir)
    listing.sort()
    for infile in listing:
        if infile[-3:] == '.gz':
            curs_list.execute( 'INSERT INTO files (file, processed) VALUES (?,?)', (infile, 0))
    conn_list.commit()
    print ("Pubchem_Listing.db created")

    # make db files
    if halogens:
        conn_pubchem = sqlite3.connect(dbs_dir + "/Pubchem_MAGMa_halo.db")
    else:
        conn_pubchem = sqlite3.connect(dbs_dir + "/Pubchem_MAGMa.db")
    curs_pubchem = conn_pubchem.cursor()
    try:
        curs_pubchem.execute(
            "CREATE TABLE molecules (cid INTEGER PRIMARY KEY, mim INTEGER NOT NULL, charge INTEGER NOT NULL, natoms INTEGER NOT NULL, molblock TEXT, inchikey TEXT, smiles TEXT, molform TEXT, name TEXT, refscore INTEGER, logp INT)")
        conn_pubchem.commit()
        print ("Pubchem_MAGMa.db created")
    except:
        print ("Pubchem_MAGMa.db already exists (or error creating it)")

    if len(kegg_info) > 0:
        if halogens:
            conn_kegg = sqlite3.connect(
                dbs_dir + "/Pubchem_MAGMa_Kegg_halo.db")
        else:
            conn_kegg = sqlite3.connect(dbs_dir + "/Pubchem_MAGMa_Kegg.db")
        curs_kegg = conn_kegg.cursor()
        try:
            curs_kegg.execute(
                "CREATE TABLE molecules (cid INTEGER PRIMARY KEY, mim INTEGER NOT NULL, charge INTEGER NOT NULL, natoms INTEGER NOT NULL, molblock TEXT, inchikey TEXT, smiles TEXT, molform TEXT, name TEXT, reference TEXT, logp INT)")
            conn_kegg.commit()
            print ("Pubchem_MAGMa_kegg.db created")
        except:
            print (
                "Pubchem_MAGMa_kegg.db already exists (or error creating it)")

    # fill databases
    ready = False
    curr_id = 0
    curr_id_kegg = 0
    memstore = {}
    memstore_kegg = {}
    while not ready:
        curs_list.execute('SELECT file FROM files WHERE processed == 0')
        try:
            filename = curs_list.fetchone()[0]
        except:
            ready = True
            continue
        starttime = time.time()
        print filename
        sdfile = gzip.open(pubchem_dir + filename)
        line = '$$$$'
        while line != "":
            record = []
            hatoms = []
            hbonds = 0
            skip = False
            halo = False
            ionized = 0
            # read heading:
            for x in range(4):
                line = sdfile.readline()
                record.append(line)
            if line == "":
                continue
            natoms = int(line[:3])
            nbonds = int(line[3:6])
            for x in range(natoms):
                line = sdfile.readline()
                if line[31:33] == 'H ':
                    # remove hydrogens
                    hatoms.append(x + 1)
                    continue
                if line[31:33] != 'H ' and len(hatoms) > 0:
                    exit('Error: heavy atoms after hydrogens' + record[0][:-1])
                if line[31:33] not in ['C ', 'N ', 'O ', 'P ', 'S ', 'F ', 'Cl', 'Br', 'I ']:
                    # filter non-organic compounds
                    skip = True
                if line[31:33] in ['F ', 'Cl', 'Br', 'I ']:
                    halo = True
                elif line[50:51] != '0':
                    # this flag has something to do with polymeric structures
                    # and resulted in deviation between calculated and given inchikeys, skip
                    skip = True
                elif line[38:39] == '4':  # radical
                    # radical, resulted in deviation between calculated and given inchikeys
                    skip = True
                record.append(line[:42] + '\n')
            for x in range(nbonds):
                line = sdfile.readline()
                if int(line[:3]) in hatoms or int(line[3:6]) in hatoms:
                    # remove bonds involving hydrogens
                    hbonds += 1
                    continue
                # use bonds with stereoflags set to zero
                record.append(line[:9] + '  0\n')
            while line != 'M  END\n' and line != '':
                line = sdfile.readline()
                record.append(line)
                if line[:6] == 'M  ISO':
                    skip = True
                    print 'Skipped isotopically labeled:', record[0][:-1]
            while line != "$$$$\n" and line != "":
                line = sdfile.readline()
                if line == "> <PUBCHEM_MONOISOTOPIC_WEIGHT>\n":
                    mim = float(sdfile.readline()[:-1])
                elif line == "> <PUBCHEM_COMPONENT_COUNT>\n":
                    comp_count = int(sdfile.readline()[:-1])
                elif line == "> <PUBCHEM_IUPAC_INCHIKEY>\n":
                    line = sdfile.readline()
                    inchikey = line[:14]
                elif line == "> <PUBCHEM_OPENEYE_CAN_SMILES>\n":
                    smiles = sdfile.readline()[:-1]
                elif line == "> <PUBCHEM_MOLECULAR_FORMULA>\n":
                    molform = sdfile.readline()[:-1]
                elif line == "> <PUBCHEM_TOTAL_CHARGE>\n":
                    charge = int(sdfile.readline()[:-1])
                elif line == "> <PUBCHEM_IUPAC_NAME>\n":
                    iupac_name = sdfile.readline()[:-1]
                elif line == "> <PUBCHEM_HEAVY_ATOM_COUNT>\n":
                    heavy_atoms = int(sdfile.readline()[:-1])
                elif line[:17] == "> <PUBCHEM_XLOGP3":
                    logp = float(sdfile.readline()[:-1])
            if line != "" and \
                    not skip and \
                    halo == halogens and \
                    mim <= 1200.0 and \
                    comp_count == 1:
                charge = 0
                if '-' in molform:
                    # exclude compounds with charge < -1
                    if molform[-1] == '-':
                        charge = -1
                    else:
                        continue
                elif '+' in molform:
                    # exclude compounds with charge > 1
                    if molform[-1] == '+':
                        charge = 1
                    else:
                        continue
                for x in ['C(=O)[O-]', 'S(=O)(=O)[O-]', '[NH+]', '[NH2+]', '[NH3+]', '[NH4+]']:
                    if smiles.find(x) >= 0:
                        ionized = 1
                record[3] = repr(natoms - len(hatoms)).rjust(3) + repr(nbonds - hbonds).rjust(3) + record[3][6:]
                cid = int(record[0][:-1])
                try:
                    curs_names.execute(
                        'SELECT name,refscore FROM names WHERE cid == ?', (cid,))
                    result = curs_names.fetchone()
                    molname = result[0]
                    refscore = int(result[1])
                except:
                    molname = iupac_name
                    refscore = 1
                molblock = base64.encodestring(zlib.compress(''.join(record)))
                if inchikey in memstore:
                    dbcid, dbrefscore, dbionized = memstore[inchikey]
                    # prefer CID's with higher refscore, then prefer non-ionized CID's
                    if (refscore > dbrefscore) or \
                            (refscore == dbrefscore and dbionized > ionized):
                        curs_pubchem.execute('''UPDATE molecules SET cid=?, mim=?, charge=?,  molblock=?, smiles=?,
                                                molform=?, name=?, refscore=?, logp=? WHERE cid == ?''', (
                                                cid,
                                                int(mim * 1e6),
                                                charge,
                                                unicode(molblock),
                                                unicode(smiles),
                                                unicode(molform),
                                                unicode(molname),
                                                refscore,
                                                int(logp * 10),
                                                dbcid))
                        memstore[inchikey] = (cid, refscore, ionized)
                else:
                    curr_id += 1
                    curs_pubchem.execute('''INSERT INTO molecules (cid,mim,charge,natoms,molblock,inchikey,smiles,
                                            molform,name,refscore,logp) VALUES (?,?,?,?,?,?,?,?,?,?,?)''', (
                                            cid,
                                            int(mim * 1e6),
                                            charge,
                                            heavy_atoms,
                                            unicode(molblock),
                                            unicode(inchikey),
                                            unicode(smiles),
                                            unicode(molform),
                                            unicode(molname),
                                            refscore,
                                            int(logp * 10)))
                    memstore[inchikey] = (cid, refscore, ionized)
                if int(cid) in kegg_info:
                    molname, keggid = kegg_info[cid]
                    if inchikey in memstore_kegg:
                        dbcid, reference, dbionized = memstore_kegg[inchikey]
                        reference = reference + ',' + keggid
                        print 'Duplicates:', reference, molname
                        # prefer non-ionized CID's
                        if dbionized > ionized:
                            curs_kegg.execute('''UPDATE molecules SET cid=?, mim=?, charge=?, molblock=?, smiles=?,
                                                 molform=?, name=?, reference=?, logp=? WHERE cid == ?''', (
                                                    cid,
                                                    int(mim * 1e6),
                                                    charge,
                                                    unicode(molblock),
                                                    unicode(molform),
                                                    unicode(smiles),
                                                    unicode(molname),
                                                    unicode(reference),
                                                    int(logp * 10),
                                                    dbcid))
                            memstore_kegg[inchikey] = (cid, reference, ionized)
                        else:
                            curs_kegg.execute('UPDATE molecules SET reference=? WHERE cid == ?', (
                                unicode(reference),
                                dbcid)
                            )
                            memstore_kegg[inchikey] = (
                                dbcid, reference, dbionized)
                    else:
                        curr_id += 1
                        curs_kegg.execute('''INSERT INTO molecules (cid, mim, charge, natoms, molblock, inchikey,
                                             smiles, molform, name, reference, logp) VALUES (?,?,?,?,?,?,?,?,?,?,?)''', (
                                                cid,
                                                int(mim * 1e6),
                                                charge,
                                                heavy_atoms,
                                                unicode(molblock),
                                                unicode(inchikey),
                                                unicode(smiles),
                                                unicode(molform),
                                                unicode(molname),
                                                unicode(keggid),
                                                int(logp * 10)))
                        memstore_kegg[inchikey] = (cid, keggid, ionized)

        conn_pubchem.commit()
        if len(kegg_info) > 0:
            conn_kegg.commit()
        elapsed_time = (time.time() - starttime) / 60
        curs_list.execute(
            'UPDATE files SET processed = 1, time = ? WHERE file = ?', (elapsed_time, filename))
        conn_list.commit()
        print filename + ": ", elapsed_time

    memstore = {}  # free up some memory
    memstore_kegg = {}
    print "Creating index ..."
    #curs_pubchem.execute('PRAGMA temp_store = 2')
    curs_pubchem.execute(
        'CREATE INDEX idx_cover ON molecules (charge,mim)')
    conn_pubchem.commit()
    if len(kegg_info) > 0:
        #curs_kegg.execute('PRAGMA temp_store = 2')
        curs_kegg.execute(
            'CREATE INDEX idx_cover ON molecules (charge,mim)')
        conn_kegg.commit()

# main
mainparser = argparse.ArgumentParser(description=__doc__)
mainparser.add_argument(
    '--version', action='version', version='%(prog)s ' + version())
subparsers = mainparser.add_subparsers(title='Sub-commands')

sc = subparsers.add_parser(
    "update", help=update.__doc__, description=update.__doc__)
sc.add_argument('-d', '--data_dir',
                help="Directory where PubChem data is stored (default: %(default)s)", default="./", type=str)
sc.set_defaults(func=update)

sc = subparsers.add_parser(
    "process", help=process.__doc__, description=process.__doc__)
sc.add_argument('-k', '--kegg', default=None, type=str, help="""File with information of Kegg compounds in PubChem.
                        Must be downloaded from http://pubchem.ncbi.nlm.nih.gov/ as follows:
                        Substance => Advanced Search => Source=KEGG => Search;
                        Display Settings: Format = ID Map => Apply;
                        Send to: File
                        (default: %(default)s)""")
sc.add_argument('-s', '--skip_names',
                help="Skip update of PubChem names db (default: %(default)s)", action="store_true")
sc.add_argument('-f', '--halogens',
                help="Generate database with halogenated compounds (default: %(default)s)", action="store_true")
sc.add_argument('-d', '--data_dir',
                help="Directory where PubChem data has been downloaded (default: %(default)s)", default="./", type=str)
sc.add_argument('-b', '--database_dir',
                help="Directory where MAGMa databases will be stored (default: %(default)s)", default="./", type=str)
sc.set_defaults(func=process)

args = mainparser.parse_args(sys.argv[1:])
args.func(args)

