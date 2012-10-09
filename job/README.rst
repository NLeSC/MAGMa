MAGMa
=====

MAGMa is a abbreviation for 'Ms Annotation based on in silico Generated Metabolites'.
Magma subproject which performs actual calculations.

Requirements
------------

MAGMa requires CDK java package via communicates with CDK using JPype python package.

Running on cluster
------------------

On the compute node not all dependencies of Magma will be installed.
By freezing the magma application on the head node we include all dependencies like rdkit.

On head node ::

   pip install bbfreeze
   python setup.py bdist_bbfreeze
   cd dist
   chmod +x dist/Magma-<version>/Magma-<version>-py2.7.egg/magma/script/reactor
   tar -zcf Magma-<version>.tar.gz Magma-<version>

On compute node ::

   tar -zxf Magma-<version>.tar.gz
   ./Magma-<version>/magma ...

