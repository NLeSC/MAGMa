#!/usr/bin/env python
"Update local HMDB database"

import os
import zlib
import sqlite3
import sys,argparse,urllib2,zipfile,StringIO
from rdkit import Chem,Geometry
from rdkit.Chem import AllChem, Descriptors

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
Hmass=mims[1]

def version():
    return '1.0'

def process_hmdb(args):
    conn = sqlite3.connect(args.database_dir+'/HMDB_MAGMa.db')
    c = conn.cursor()
    try:
        c.execute("CREATE TABLE molecules (id TEXT PRIMARY KEY, mim INTEGER NOT NULL, charge INTEGER NOT NULL, natoms INTEGER NOT NULL, molblock BLOB, inchikey TEXT, smiles TEXT, molform TEXT, name TEXT, reference TEXT, logp INT)")
        conn.commit()
        print ("HMDB_MAGMa.db created")
    except:
        print ("HMDB_MAGMa.db already exists (or error creating it)")
        exit()

    if args.data_dir == None:
        zf=urllib2.urlopen('http://www.hmdb.ca/downloads/structures.zip')
    else:
        zf=open(args.data_dir+'structures.zip')
    sdfile=zipfile.ZipFile(StringIO.StringIO(zf.read())).open('structures.sdf')
    
    memstore={}
    line='$$$$'
    while line!="":
        record=[]
        amap={}
        skip=False
        ionized=0
        #read heading:
        for x in range(4):
            line=sdfile.readline()
            #print line[:-1]
            record.append(line)
        if line == "":
            continue
        natoms=int(record[-1][:3])
        nbonds=int(record[-1][3:6])
        bonds=0
        y=0
        for x in range(natoms):
            line=sdfile.readline()
            #print line[:-1]
            if line[31:33] == 'H ':
                #print x,":",line,
                continue # remove hydrogens
            y+=1
            amap[x+1]=y
            #if line[31:33] != 'H ' and len(hatoms)>0:
            #    print 'Skipped: heavy atoms after hydrogens --- ' # +str(record)
            #    skip=True
            #    # exit()
            if line[31:33] not in ['C ','N ','O ','P ','S ','F ','Cl','Br','I ']:
                skip=True # filter non-organic compounds
            elif line[50:51] != '0': # this flag seems to have something to do with polymeric structures
                skip=True # -> resulted in deviation between calculated and given inchikeys
            elif line[38:39] == '4': # radical
                skip=True # -> resulted in deviation between calculated and given inchikeys
            record.append(line[:42]+'\n')
        #print amap
        for x in range(nbonds):
            line=sdfile.readline()
            #print line[:-1]
            a1=int(line[:3])
            a2=int(line[3:6])
            if a1 in amap and a2 in amap: # skip bonds involving hydrogens
                bonds+=1
                #record.append(line[:9]+'  0\n') # use bonds with stereoflags set to zero
                record.append('%3i%3i%s  0\n' % (amap[a1],amap[a2],line[6:9])) # use bonds with stereoflags set to zero
        while line != 'M  END\n' and line != '':
            line=sdfile.readline()
            record.append(line)
            if line[:6]=='M  ISO':
                 skip=True
                 print 'Skipped isotopically labeled:',record[0][:-1]
        while line != "$$$$\n" and line != "":
            line=sdfile.readline()
            if line == "> <HMDB_ID>\n":
                hmdb_id=str(sdfile.readline()[:-1])
            if line == "> <GENERIC_NAME>\n":
                molname=str(sdfile.readline()[:-1])
        if line != "" and skip==False:
            record[3]=repr(y).rjust(3)+repr(bonds).rjust(3)+record[3][6:]
            molblock=''.join(record)
            mol=Chem.MolFromMolBlock(molblock)
            if mol==None or mol.GetNumAtoms()==0:
                continue
            smiles=Chem.MolToSmiles(mol)
            if len(Chem.GetMolFrags(mol))>1:
                print 'complex:',hmdb_id,smiles
                continue
            mass=0.0
            for atom in mol.GetAtoms():
                mass+=mims[atom.GetAtomicNum()]+Hmass*(atom.GetNumImplicitHs()+atom.GetNumExplicitHs())
            #if skip:
            #    print 'non-organic element in:',hmdb_id,smiles
            #    continue
            #print smiles, hmdb_id
            inchi=Chem.AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))
            conf = mol.GetConformer(0)
            molblock = zlib.compress(''.join(record))
            molform=Chem.rdMolDescriptors.CalcMolFormula(mol)
            mim=Chem.rdMolDescriptors.CalcExactMolWt(mol)
            charge=0
            if '-' in molform:
                if molform[-1]=='-':
                    charge=-1
                    # mim+=0.0005486 # account for mass of electron
                else:
                    continue
            elif '+' in molform:
                if molform[-1]=='+':
                    charge=1
                    # mim-=0.0005486 # account for mass of electron
                else:
                    continue
            if mim > 1200.0:
                print 'molecule to heavy:',hmdb_id,smiles
                continue
            natoms=mol.GetNumHeavyAtoms()
            logp=Chem.Crippen.MolLogP(mol)
            inchikey = Chem.AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))[:14]
            ionized=0
            for x in ['C(=O)[O-]','[NH+]','[NH2+]','[NH3+]','[NH4+]']:
                if smiles.find(x)>=0:
                    ionized=1    
            if inchikey in memstore:
                dbid,reference,dbionized=memstore[inchikey]
                reference=reference+','+hmdb_id
                print 'Duplicates:',reference,molname
                if dbionized > ionized: # prefer non-ionized CID's
                    c.execute('UPDATE molecules SET id=?, mim=?, charge=?, molblock=?, smiles=?, molform=?, name=?, reference=?, logp=? WHERE id == ?', (
                                                             hmdb_id,
                                                             int(mim*1e6),
                                                             charge,
                                                             buffer(molblock),
                                                             unicode(smiles),
                                                             unicode(molform),
                                                             unicode(molname, 'utf-8', 'xmlcharrefreplace'),
                                                             unicode(reference),
                                                             int(logp*10),
                                                             dbid)
                                                             )
                    memstore[inchikey]=(hmdb_id,reference,ionized)
                else:
                    c.execute('UPDATE molecules SET reference=? WHERE id == ?', (
                                                             unicode(reference),
                                                             dbid)
                                                             )
                    memstore[inchikey]=(dbid,reference,dbionized)
            else:
                c.execute('INSERT INTO molecules (id,mim,charge,natoms,molblock,inchikey,smiles,molform,name,reference,logp) VALUES (?,?,?,?,?,?,?,?,?,?,?)', (
                                                             hmdb_id,
                                                             int(mim*1e6),
                                                             charge,
                                                             int(natoms),
                                                             buffer(molblock),
                                                             unicode(inchikey),
                                                             unicode(smiles),
                                                             unicode(molform),
                                                             unicode(molname, 'utf-8', 'xmlcharrefreplace'),
                                                             unicode(hmdb_id),
                                                             int(logp*10)
                                                             )
                      )
                memstore[inchikey]=(hmdb_id,hmdb_id,ionized)
    conn.commit()

    print "Creating index ..."
    c.execute('PRAGMA temp_store = 2')
    c.execute('CREATE INDEX idx_cover ON molecules (charge,mim,natoms,reference,molform,inchikey,smiles,name,molblock,logp)')
    conn.commit()

#main
mainparser = argparse.ArgumentParser(description=__doc__)
mainparser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version())
mainparser.add_argument('-d', '--data_dir', help="Directory where HMDB structures.zip file is stored (default: %(default)s, attempt to read directly from HMDB server)", default=None,type=str)
mainparser.add_argument('-b', '--database_dir', help="Directory where HMDB database is stored (default: %(default)s)", default="./",type=str)
mainparser.set_defaults(func=process_hmdb)

args = mainparser.parse_args(sys.argv[1:])
args.func(args)
