Create local PubChem and HMDB structure databases
=================================================

PubChem is a large database of chemical molecules and their activities against biological assays.
HMDB is a database of human metabolites
For more information see http://pubchem.ncbi.nlm.nih.gov/ and http://www.hmdb.ca/

Create local PubChem and Kegg (=subset of PubChem) databases:
-------------------------------------------------------------

1) Download/update local copy of PubChem data files:

.. code-block::

   usage: process_pubchem.py update [-h] [-d DATA_DIR]
 
   Update downloads from PubChem server
   
   optional arguments:
     -h, --help            show this help message and exit
     -d DATA_DIR, --data_dir DATA_DIR
                           Directory where PubChem data is stored (default: ./)

2) Manual download of file with information of Kegg compounds in PubChem.
   Must be downloaded from http://pubchem.ncbi.nlm.nih.gov/ as follows:
   Substance => Advanced => SourceName: kegg => Search;
   Display Settings: Format = ID Map => Apply;
   Send to: File

3) Process the PubChem data twice, once for non-halogenated compounds and
   once for halogenated compounds (-f option)

.. code-block::

   usage: process_pubchem.py process [-h] [-k KEGG] [-s] [-f] [-d DATA_DIR]
                                     [-b DATABASE_DIR]
   
   Update local PubChem databases. This has to be run twice, with and without -f option
   to generate db's for halogenated and non-halogenated compounds respectively
   
   optional arguments:
     -h, --help            show this help message and exit
     -k KEGG, --kegg KEGG  File obtained in step 2
     -s, --skip_names      Skip update of PubChem names db (default: False)
     -f, --halogens        Generate database with halogenated compounds (default:
                           False)
     -d DATA_DIR, --data_dir DATA_DIR
                           Directory where PubChem data is stored (default: ./)
     -b DATABASE_DIR, --database_dir DATABASE_DIR
                           Directory where PubChem databases are stored (default:
                           ./)

Create local HMDB database:
---------------------------

.. code-block::

   usage: process_hmdb.py [-h] [-v] [-d DATA_DIR] [-b DATABASE_DIR]

   Update local HMDB database
   
   optional arguments:
     -h, --help            show this help message and exit
     -v, --version         show program's version number and exit
     -d DATA_DIR, --data_dir DATA_DIR
                           Directory where HMDB structures.zip file is stored
                           (default: None, attempt to read directly from HMDB server)
     -b DATABASE_DIR, --database_dir DATABASE_DIR
                           Directory where HMDB database is stored (default: ./)
   