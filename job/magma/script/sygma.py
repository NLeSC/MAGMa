import fileinput,sys,getopt,base64,subprocess,codecs
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Geometry
from rdkit.Chem import Descriptors
import magma

def usage():
   sys.stderr.write("""
   Usage: metabolize -t type -i input [-n nsteps] [-m metabolism] [-h] output.sqlite
   
      -o: output SDF
      -i: input = smilesstring, filename or base64-encoded ctab
      -m: metabolism = 1 and/or 2 for phase 1 and 2 biotransformations
      -n: nsteps
      -h: print this help
   """+"\n")

def main(argv=sys.argv):
   # Set Defaults
   nsteps=2
   metabolism="12"
   # Process input parameters
   try:
      opts,args = getopt.getopt(sys.argv[1:], "o:i:n:m:h")
   except getopt.GetoptError, err:
      sys.stderr.write(str(err))
      usage()
      sys.exit(2)

   parentmol=[]
   for o,a in opts:
      if o=="-h":
         usage()
         sys.exit(2)
      elif o=="-i":
         inputline=a.split("|")
         if inputline[0]=="smiles":
               parentmol.append(Chem.MolFromSmiles(inputline[1]))
               AllChem.Compute2DCoords(parentmol[-1])
               parentmol[-1].SetProp("_Name",inputline[2])
         elif inputline[0]=="smilesfile":
            smilesfile=codecs.open(inputline[1],'r',encoding='utf-8')
            for line in smilesfile:
               if line[:-1] != "":
                  print line[:-1].encode('ascii','xmlcharrefreplace')
                  parentmol.append(Chem.MolFromSmiles(line.split("|")[0].encode('ascii')))
                  AllChem.Compute2DCoords(parentmol[-1])
                  parentmol[-1].SetProp("_Name",line[:-1].split("|")[1].encode('ascii','xmlcharrefreplace'))
            smilesfile.close()
         elif inputline[0]=="molfile":
            parentmol.append(Chem.MolFromMolFile(inputline[1]))
         elif inputline[0]=="molstruct":
            parentmol.append(Chem.MolFromMolBlock(base64.standard_b64decode(inputline[1])))
         inputmol=a
      elif o=="-n":
         nsteps=int(a)
      elif o=="-m":
         metabolism=a
      elif o=="-o":
         output=a

   print len(parentmol)

   if len(args) > 0:
      magma.set_DB(args[0])
   else:
      usage()
      sys.exit(2)

   for x in parentmol:
      if x.GetNumAtoms()>1: # Chemdoodle editor outputs methane by default
         mol=Chem.MolToMolBlock(x)
         magma.add_metabolite(mol,mol.split('\n')[0],1.0,0,'PARENT',1)

   if nsteps>0:
      magma.metabolize_all(metabolism,nsteps)
   magma.commit_DB()
   if output!="":
      magma.printSDF(open(a,'w'))
