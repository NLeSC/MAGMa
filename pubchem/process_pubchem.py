#!/usr/bin/env python

import os
import zlib,gzip
import sqlite3
import time

inputdir="/media/PubChem/ftp.ebi.ac.uk/pub/databases/pubchem/Compound/CURRENT-Full/SDF/"

conn2 = sqlite3.connect("Pubchem_Names.db")
c2 = conn2.cursor()
try:
   c2.execute("CREATE TABLE names (cid INTEGER PRIMARY KEY, name TEXT, refscore INTEGER)")
   #currently using a Synonym db created with Stefans java script
   #could be generated with something like:
   namesfile=gzip.open('CID-Synonym-filtered.gz')
   sidfile=gzip.open('CID-SID.gz')
   splitsid=sidfile.readline().split("\t")
   splitnames=namesfile.readline().split("\t")
   while splitsid != ['']:
      curr_cid=splitsid[0]
      sid_count=1
      names_count=1
      name=''
      splitsid=sidfile.readline().split("\t")
      while splitsid[0]==curr_cid:
         sid_count+=1
         splitsid=sidfile.readline().split("\t")
      if splitnames[0]==curr_cid:
         name=splitnames[1][:-1]
      while splitnames[0]==curr_cid:
         names_count+=1
         splitnames=namesfile.readline().split("\t")
      c2.execute('INSERT INTO names (cid, name, refscore) VALUES (?,?,?)', (int(curr_cid),name,sid_count*names_count))
      if names_count==1:
         print curr_cid,name,sid_count,names_count,sid_count*names_count
   conn2.commit()
except:
   print ("Pubchem_Names.db already exists (or error creating it)")

conn3 = sqlite3.connect('Pubchem_Listing.db')
c3=conn3.cursor()
try:
   c3.execute("CREATE TABLE files (file TEXT UNIQUE NOT NULL, processed INTEGER, time REAL)")
   conn3.commit()
   listing = os.listdir(inputdir)
   listing.sort()
   for infile in listing:
      if infile[-3:] == '.gz':
         c3.execute('INSERT INTO files (file, processed) VALUES (?,?)',(infile,0))
   conn3.commit()
   print ("Pubchem_Listing.db created")
except:
   print ("Pubchem_Listing.db already exists (or error creating it)")

conn = sqlite3.connect("Pubchem_MAGMa_new.db")
c = conn.cursor()
try:
   c.execute("CREATE TABLE molecules (cid INTEGER PRIMARY KEY, mim INTEGER NOT NULL, natoms INTEGER NOT NULL, molblock BLOB, inchikey TEXT, molform TEXT, name TEXT, refscore INTEGER, logp INT)")
   conn.commit()
   print ("Pubchem_MAGMa.db created")
except:
   print ("Pubchem_MAGMa.db already exists (or error creating it)")

ready=False
curr_id=0
memstore={}
while not ready:
   c3.execute('SELECT file FROM files WHERE processed == 0')
   try:
      filename=c3.fetchone()[0]
   except:
      ready=True
      continue
   starttime=time.time()
   sdfile=gzip.open(inputdir+filename)
   record=[]
   hatoms=[]
   hbonds=0
   skip=False
   ionized=0
   line='$$$$'
   while line!="":
      #read heading:
      for x in range(4):
         line=sdfile.readline()
         record.append(line)
      if line == "":
         continue
      natoms=int(record[-1][:3])
      nbonds=int(record[-1][3:6])
      for x in range(natoms):
         line=sdfile.readline()
         if line[31:33] == 'H ':
            hatoms.append(x+1)
            continue # remove hydrogens
         if line[31:33] != 'H ' and len(hatoms)>0:
            exit('Error: heavy atoms after hydrogens')
         if line[31:33] not in ['C ','N ','O ','P ','S ']:
            skip=True # filter non-organic compounds
         elif line[50:51] != '0': # this flag seems to have something to do with polymeric structures
            skip=True # -> resulted in deviation between calculated and given inchikeys
         elif line[38:39] == '4': # radical
            skip=True # -> resulted in deviation between calculated and given inchikeys
         record.append(line[:42]+'\n')
      for x in range(nbonds):
         line=sdfile.readline()
         if int(line[:3]) in hatoms or int(line[3:6]) in hatoms:
            hbonds+=1
            continue # remove bonds involving hydrogens
         record.append(line[:12]+'\n')
      while line != 'M  END\n' and line != '':
         line=sdfile.readline()
         record.append(line)
      while line != "$$$$\n" and line != "":
         line=sdfile.readline()
         if line == "> <PUBCHEM_MONOISOTOPIC_WEIGHT>\n":
            mim=float(sdfile.readline()[:-1])
         elif line == "> <PUBCHEM_COMPONENT_COUNT>\n":
            comp_count=int(sdfile.readline()[:-1])
         elif line == "> <PUBCHEM_IUPAC_INCHIKEY>\n":
            line=sdfile.readline()
            #print line[:-1],
            inchikey=line[:14]
         elif line == "> <PUBCHEM_OPENEYE_CAN_SMILES>\n":
            smiles=sdfile.readline()[:-1]
         elif line == "> <PUBCHEM_MOLECULAR_FORMULA>\n":
            molform=sdfile.readline()[:-1]
         elif line == "> <PUBCHEM_TOTAL_CHARGE>\n":
            charge=int(sdfile.readline()[:-1])
         elif line == "> <PUBCHEM_IUPAC_NAME>\n":
            iupac_name=sdfile.readline()[:-1]
         elif line == "> <PUBCHEM_HEAVY_ATOM_COUNT>\n":
            heavy_atoms=int(sdfile.readline()[:-1])
         elif line[:17] == "> <PUBCHEM_XLOGP3":
            logp=float(sdfile.readline()[:-1])
      if line != "" and not skip and \
            mim <= 1200.0 and \
            comp_count == 1:
         for x in ['C(=O)[O-]','S(=O)(=O)[O-]','[NH+]','[NH2+]','[NH3+]','[NH4+]']:
            if smiles.find(x)>=0:
               ionized=1
         record[3]=repr(natoms-len(hatoms)).rjust(3)+repr(nbonds-hbonds).rjust(3)+record[3][6:]
         cid=int(record[0][:-1])
         try:
            c2.execute('SELECT name,refscore FROM names WHERE cid == ?',(cid,))
            result=c2.fetchone()
            molname=result[0]
            refscore=int(result[1])
         except:
            molname=iupac_name
            refscore=1
         if inchikey in memstore:
            dbcid,dbrefscore,dbionized=memstore[inchikey]
            if (refscore > dbrefscore) or \
               (refscore == dbrefscore and dbionized > ionized): # prefer CID's with higher refscore, then prefer non-ionized CID's
               molblock = zlib.compress(''.join(record))
               c.execute('UPDATE molecules SET cid=?, molblock=?, name=?, refscore=?, logp=? WHERE cid == ?', (
                                                            cid,
                                                            buffer(molblock),
                                                            unicode(molname),
                                                            refscore,
                                                            int(logp*10),
                                                            dbcid)
                                                            )
               memstore[inchikey]=(cid,refscore,ionized)
         else:
            molblock = zlib.compress(''.join(record))
            curr_id+=1
            c.execute('INSERT INTO molecules (cid,mim,natoms,molblock,inchikey,molform,name,refscore,logp) VALUES (?,?,?,?,?,?,?,?,?)', (
                                                            cid,
                                                            int(mim*1e6),
                                                            heavy_atoms,
                                                            buffer(molblock),
                                                            unicode(inchikey),
                                                            unicode(molform),
                                                            unicode(molname),
                                                            refscore,
                                                            int(logp*10)
                                                            )
                     )
            memstore[inchikey]=(cid,refscore,ionized)
      record=[]
      hatoms=[]
      hbonds=0
      skip=False
      ionized=0
   conn.commit()
   elapsed_time=(time.time()-starttime)/60
   c3.execute('UPDATE files SET processed = 1, time = ? WHERE file = ?',(elapsed_time,filename))
   conn3.commit() 
   print filename+": ",elapsed_time

print "Creating index ..."
c.execute('PRAGMA temp_store = 2')
c.execute('CREATE INDEX idx_cover ON molecules (mim,natoms,refscore,molform,inchikey,name,molblock,logp)')
conn.commit()

