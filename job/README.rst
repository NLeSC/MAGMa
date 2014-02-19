MAGMa
=====

MAGMa is a abbreviation for 'Ms Annotation based on in silico Generated Metabolites'.
Magma subproject which performs actual calculations.

Requirements
------------

MAGMa requires
* RDKit with INCHI support.
* libxml2-dev and libxslt-dev for mzxml file parsing.
* CactVS. See http://www2.ccc.uni-erlangen.de/software/cactvs/

RDKit installation
~~~~~~~~~~~~~~~~~~

See http://code.google.com/p/rdkit/wiki/BuildingWithCmake for RDKit installation instructions.

.. code-block:: bash

   sudo apt-get install flex bison build-essential python-numpy cmake python-dev sqlite3 libsqlite3-dev libboost-dev libboost-python-dev libboost-regex-dev
   virtualenv env
   . env/bin/activate
   pip install numpy
   tar -zxf RDKit_2013_09_2.tgz
   cd RDKit_2013_09_2/
   cd ../External/INCHI-API/
   mkdir build
   cd build
   cmake -DRDK_BUILD_INCHI_SUPPORT=ON ..
   make -j 4
   make install

To `env/bin/activate` add:

.. code-block:: bash

   export RDBASE=<somedir>/RDKit_2013_09_2
   export LD_LIBRARY_PATH=$RDBASE/lib
   export PYTHONPATH=$PYTHONPATH:$RDBASE

`<somedir>` should be replaced with path where RDKit was untarred.

Cactvs installation
~~~~~~~~~~~~~~~~~~~

1. On http://www.xemistry.com/ follow the Academic downloads link.
2. Download the tarball for your OS.

.. code-block:: bash

   cd /tmp
   mkdir v
   cd v
   tar -zxf ../cactvstools-Linux3.1-SuSE12.1-64-3.420.tar.gz
   ./installme
   # use defaults

MAGMa requires the path where Cactvs was installed.
For example during `./installme` of Cactvs version `3.420`, filled in `/usr/local/lib` for the libraries location then MAGMa
requires `/usr/local/lib/cactvs3.420`.

If metabolism is not working compare `magma/script/csreact` with `csbr` (available in Cactvs bin dir).
The `magma/script/csreact` must use the same paths/versions as the `csbr` script.

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

Configuration is optional.
Defaults to using 'rdkit' as chemical engine.

A 'magma_job.ini' config file is read from current working directory or from users home directory (~/).

Exampe config file ::

.. code-block:: INI

   [magma job]
   # Location of structure database to fetch candidate molecules to match against ms peak trees
   # db is expected to be available at where job is executed
   structure_database.pubchem = /media/PubChem/Pubchem_MAGMa.db
   structure_database.kegg = /media/PubChem/Pubchem_MAGMa_kegg.db
   structure_database.hmdb = /home/ridderl/hmdb/HMDB_MAGMa.db

   # MACS authentication, used for sending progress reports to MAGMa web application
   macs.id = <MAC key identifier>
   macs.key = <MAC key>

   # choose cdk or rdkit as chemical engine
   chemical_engine = rdkit

   # choose cactvs or reactor as metabolism engine
   metabolism_engine = cactvs

   # Cactvs is needed to perform reactions
   cactvs_root = /usr/local/lib/cactvs3.409

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

