MAGMa
=====

MAGMa is a abbreviation for 'Ms Annotation based on in silico Generated Metabolites'.
Magma subproject which performs actual calculations.

Requirements
------------

MAGMa requires RDKit with INCHI support.
See http://code.google.com/p/rdkit/wiki/BuildingWithCmake for RDKit installation instructions.

Development installation
------------------------

.. code-block:: bash

   python setup.py develop

Usage
-----

Annotate a tree file using PubChem database:

.. code-block:: bash

   echo '353.087494: 69989984 (191.055756: 54674544 (85.029587: 2596121, 93.034615: 1720164, 109.029442: 917026, 111.045067: 1104891 (81.034691: 28070, 83.014069: 7618, 83.050339: 25471, 93.034599: 36300, 96.021790: 8453), 127.039917: 2890439 (57.034718: 16911, 81.034706: 41459, 83.050301: 35131, 85.029533: 236887, 99.045074: 73742, 109.029404: 78094), 171.029587: 905226, 173.045212: 2285841 (71.013992: 27805, 93.034569: 393710, 111.008629: 26219, 111.045029: 339595, 137.024292: 27668, 155.034653: 145773), 191.055725: 17000514), 353.087097: 4146696)' > example.tree
   magma read_ms_data --ms_data_format tree -l 5 -a 0  example.tree results.db
   magma annotate -p5 -q0 -c0 -d0 -b3 -i -1 -s pubchem -o ../pubchem/Pubchem_MAGMa_new.db,0,9999 -f results.db

Running on cluster
------------------

On the compute node not all dependencies of Magma will be installed.
By freezing the magma application on the head node we include all dependencies like rdkit.

On head node ::

.. code-block:: bash

   pip install bbfreeze
   python setup.py bdist_bbfreeze
   cd dist
   chmod +x dist/Magma-<version>/Magma-<version>-py2.7.egg/magma/script/reactor
   tar -zcf Magma-<version>.tar.gz Magma-<version>

On compute node ::

.. code-block:: bash

   tar -zxf Magma-<version>.tar.gz
   ./Magma-<version>/magma ...

