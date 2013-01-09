PubChem processing
==================

Creates a mass lookup database from PubChem SD-files.
PubChem is a database of chemical molecules and their activities against biological assays.
For more about PubChem see http://pubchem.ncbi.nlm.nih.gov/

Usage
-----

Download auxiliary files:

.. code-block:: bash

   wget ftp://ftp.ebi.ac.uk/pub/databases/pubchem/Compound/Extras/CID-Synonym-filtered.gz
   wget ftp://ftp.ebi.ac.uk/pub/databases/pubchem/Compound/Extras/CID-SID.gz

Download SD-files (approx. 36Gb):

.. code-block:: bash

   wget -m ftp://ftp.ebi.ac.uk/pub/databases/pubchem/Compound/CURRENT-Full/SDF

Convert SD-files to molecules db (Pubchem_MAGMa_new.db)
and auxiliary files to helper databases (Pubchem_Names.db, Pubchem_Listing.db):

.. code-block:: bash

   python process_pubchem.py ftp.ebi.ac.uk/pub/databases/pubchem/Compound/CURRENT-Full/SDF/
