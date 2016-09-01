The eMetabolomics research project
==================================

.. image:: https://travis-ci.org/NLeSC/MAGMa.svg?branch=master
    :target: https://travis-ci.org/NLeSC/MAGMa

.. image:: https://landscape.io/github/NLeSC/MAGMa/master/landscape.svg?style=flat
    :target: https://landscape.io/github/NLeSC/MAGMa/master
    :alt: Code Health

.. image:: https://coveralls.io/repos/NLeSC/MAGMa/badge.svg?branch=master
    :target: https://coveralls.io/r/NLeSC/MAGMa?branch=master

.. image:: https://img.shields.io/badge/docker-ready-blue.svg
    :target: https://hub.docker.com/r/nlesc/magma

The eMetabolomics project is funded by the Netherlands eScience Center and is carried out at Wageningen University and the Netherlands eScience Center in collaboration with the Netherlands Metabolomics Centre. The project develops chemo-informatics based methods for metabolite identification and biochemical network reconstruction in an integrative metabolomics data analysis workflow.

Homepage: http://www.emetabolomics.org

MAGMa is a abbreviation for 'Ms Annotation based on in silico Generated Metabolites'.

  .. image:: web/magmaweb/static/img/metabolites.png
     :alt: Screenshot MAGMa results page

Subprojects:

- emetabolomics_site - The http://www.emetabolomics.org website
- job - Runs MAGMa calculation
- joblauncher - Webservice to execute jobs
- pubchem - Processing of PubChem database, used to find mass candidates
- web - Web application to start jobs and view results

Subproject interdependencies
----------------------------

- The `emetabolomics_site` website can be used as starting pages for the `web` application.
- The `job` calculation requires a pubchem lookup database which can be made using the `pubchem` application.
- The `web` application starts `job` calculations via the `joblauncher` webservice.

Joblauncher submodule
---------------------

Use following command to initialize and fetch the joblauncher submodule:

.. code-block:: bash

    git submodule update --init

License
-------

MAGMa is released under the Apache License Version 2.0.
The MAGMa web application uses ExtJS GPLv3 with application exception.
