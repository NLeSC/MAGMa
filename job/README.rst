MAGMa
=====

MAGMa is a abbreviation for 'Ms Annotation based on in silico Generated Metabolites'.
MAGMa subproject which performs actual calculations.

Docker
------

The MAGMa command-line script can be executed in a Docker container. For more information on installing and running Docker see https://www.docker.com/.
With the Docker image it is particularly straightforward to run the "light" subcommand which takes a single MS/MS spectrum as input, retrieves candidate structures from the online databases or a local file and writes the result as smiles or sdf to stdout.

Examples:

.. code-block:: bash

   $ docker run --rm nlesc/magma light -h
   $ cat >glutathione.mgf 
   BEGIN IONS
   TITLE=CASMI 2014, Challenge 9
   PEPMASS=308.0912 100.0
   116.0165 3.2
   144.0114 6.3
   162.0219 40.2
   179.0485 100.0
   233.0590 21.6
   290.0802 5.1
   END IONS
   ^d
   $ docker run --rm -v $PWD:/data nlesc/magma light -f mgf -s hmdb glutathione.mgf
   


   

Requirements
------------

MAGMa requires
* RDKit with INCHI support.
* libxml2-dev and libxslt-dev for mzxml file parsing.

RDKit installation
~~~~~~~~~~~~~~~~~~

See https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md for RDKit installation instructions.
Download release from https://github.com/rdkit/rdkit/releases

.. code-block:: bash

   sudo apt-get install flex bison build-essential python-numpy cmake python-dev sqlite3 libsqlite3-dev libboost-dev libboost-python-dev libboost-regex-dev
   virtualenv env
   . env/bin/activate
   pip install numpy
   tar -zxf Release_2014_09_2.tar.gz
   cd rdkit-Release_2014_09_2/
   cd External/INCHI-API/
   ./download-inchi.sh
   cd ../..
   mkdir build
   cd build
   cmake -DRDK_BUILD_INCHI_SUPPORT=ON ..
   make -j 4
   make install

To `env/bin/activate` add:

.. code-block:: bash

   export RDBASE=<somedir>/rdkit-Release_2014_09_2
   export LD_LIBRARY_PATH=$RDBASE/lib
   export PYTHONPATH=$PYTHONPATH:$RDBASE

`<somedir>` should be replaced with path where RDKit was untarred.

Development installation
------------------------

.. code-block:: bash

   pip install cython
   python setup.py develop

Usage
-----

Annotate a tree file using PubChem database:

.. code-block:: bash

   echo '353.087494: 69989984 (191.055756: 54674544 (85.029587: 2596121, 93.034615: 1720164, 109.029442: 917026, 111.045067: 1104891 (81.034691: 28070, 83.014069: 7618, 83.050339: 25471, 93.034599: 36300, 96.021790: 8453), 127.039917: 2890439 (57.034718: 16911, 81.034706: 41459, 83.050301: 35131, 85.029533: 236887, 99.045074: 73742, 109.029404: 78094), 171.029587: 905226, 173.045212: 2285841 (71.013992: 27805, 93.034569: 393710, 111.008629: 26219, 111.045029: 339595, 137.024292: 27668, 155.034653: 145773), 191.055725: 17000514), 353.087097: 4146696)' > example.tree
   magma read_ms_data --ms_data_format tree -l 5 -a 0  example.tree results.db
   magma annotate -p5 -q0 -c0 -d0 -b3 -i -1 -s pubchem -o ../pubchem/Pubchem_MAGMa_new.db,0,9999 -f results.db

Configuration
-------------

A 'magma_job.ini' config file is read from users home directory (~/).

Exampe config file to read candidate molecules from the emetabolomics server:

.. code-block:: INI

   [magma job]
   # Retrieve candidate molecules from the emetabolomics server
   structure_database.online = True
   structure_database.service = http://www.emetabolomics.org/magma/molecules

Example config file to read candidate molecules from local databases (can be created by the scripts in MAGMa/pubchem):

.. code-block:: INI

   [magma job]
   # Location of structure database from which to retrieve candidate molecules locally
   structure_database.online = False
   structure_database.pubchem = /home/user/magma_databases/Pubchem_MAGMa.db
   structure_database.pubchem_halo = /home/user/magma_databases/Pubchem_MAGMa_halo.db
   structure_database.kegg = /home/user/magma_databases/Pubchem_MAGMa_kegg.db
   structure_database.kegg_halo = /home/user/magma_databases/Pubchem_MAGMa_kegg_halo.db
   structure_database.hmdb = /home/user/magma_databases/HMDB_MAGMa.db

   # MACS authentication, used for sending progress reports to MAGMa web application
   macs.id = <MAC key identifier>
   macs.key = <MAC key>

Running on cluster
------------------

On the compute node not all dependencies of Magma will be installed.
By freezing the magma application on the head node we include all dependencies like rdkit.

On head node:

.. code-block:: bash

   pip install bbfreeze
   python setup.py bdist_bbfreeze
   cd dist
   chmod +x dist/Magma-<version>/Magma-<version>-py2.7.egg/magma/script/reactor
   tar -zcf Magma-<version>.tar.gz Magma-<version>

On compute node:

.. code-block:: bash

   tar -zxf Magma-<version>.tar.gz
   ./Magma-<version>/magma ...
