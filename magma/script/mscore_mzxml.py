import sys,getopt
from rdkit import Chem
import magma

"""
RDkit dependencies:
determine bondtype
calculate MIM of parent
count implicit/explicit hydrogens (to calculate MIM of fragments)

old timing: user 18m38
"""

# import zlib

precision=0.001   # precision to identify masses
Nbonds = 4       # Allowed number of bond breaks
MSfilter = 2e5   # Intensity cutoff (absolute) to select ions
MSMSfilter = 0.1 # Intensity cutoff (relative to basePeakIntensity) to select fragment peaks
# MSfilter = 1e6   # Intensity cutoff (absolute) to select ions
# MSMSfilter = 0.5 # Intensity cutoff (relative to basePeakIntensity) to select fragment peaks
fpp = 1e-7       # precision for matching floating point values
useFragments=True # Assign fragment data
ionisation = 1   # positive ionisation mode
useMSMSonly=True

def usage():
   sys.stderr.write("""
   Usage: mscore_mzxml -f mzXMLfile [-e] [-h] [-p precision] [SyGMa_SDfile]

      -f: mzXMLfile
      -h: print this help text
      -p: precision in Dalton, default = 0.01
      -c: cutoff value to filter MS peaks (absolute)
      -d: cutoff value to filter MSMS peaks (relative to basepeak
      -d: database file (sqlite)
      -i: ionisation mode = "neg" or "pos"
      -s: skip interpretation of fragment spectra
   """+"\n")


def main(argv=sys.argv):
        sys.stderr.write('Running MScore_mzXML.py\n')

        # Process input parameters
        try:
           opts,args = getopt.getopt(argv[1:], "sef:hp:c:d:i:")
        except getopt.GetoptError, err:
           sys.stderr.write(str(err))
           usage()
           sys.exit(2)

        for o,a in opts:
           if o=="-h":
              usage()
              sys.exit(2)
           elif o=="-s":
              useFragments=False
           elif o=="-f":
              mzxmlfile=a
           elif o=="-p":
              precision=float(a)
           elif o=="-c":
              MSfilter=float(a)
           elif o=="-d":
              database=a
           elif o=="-i":
              if a == "neg":
                 ionisation=-1
              elif a == "pos":
                 ionisation=1
              else:
                 usage()
                 sys.exit(2)

        if len(args) > 0:
           magma.set_DB(args[0])
        else:
           usage()
           sys.exit(2)

        # Put run parameters in database
        magma.add_run_data(0,0,0,mzxmlfile,ionisation,useFragments,MSfilter,MSMSfilter,precision,useMSMSonly)

        # Read mzXML file and store data in database
        magma.storeMZxmlFile(mzxmlfile)
        magma.commit_DB()

        # Read spectral data from database
        magma.buildspectra()

        # Process molecules from database and store scores and fragments
        magma.searchAllMetabolites()

        magma.commit_DB()


