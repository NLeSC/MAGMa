====================
Web service consumer
====================

MAGMa has a nice web user interface.
The user interface recieves/sends data using a JSON web service.
The JSON web service can also be used by other applications.
Below is an example how to use the web service using `curl <http://curl.haxx.se/>`_, a command line HTTP client.

.. contents::

Start session
=============

.. code-block:: bash

   curl -L -c cookie.jar -b cookie.jar http://www.emetabolomics.org/magma

Returns the html for the MAGMa home page.
The output can be ignored. This was just to get a session cookie into the cookie.jar.

Submit job
==========

.. code-block:: bash

   curl -L -c cookie.jar -b cookie.jar \
   -F ms_data_file=@input.mtree -F ms_data_format=mass_tree \
   -F structure_database=pubchem \
   -F max_mz=1200 -F min_refscore=1 \
   -F ionisation_mode=-1 \
   -F max_broken_bonds=3 -F max_water_losses=1 \
   -F mz_precision=5 \
   http://www.emetabolomics.org/magma/start

Parameters:

- ``ms_data_file``, file upload of mass spectra data
- ``ms_data_format``, which format ``ms_data_file`` is in. See http://www.emetabolomics.org/magma/help for examples.
- ``structure_database``, Retrieve molecules from a database. Can be pubchem, kegg, hmdb.
- ``max_mz``, Mass filter for structure database
- ``min_refscore``, PubChem reference score filter for pubchem structure database.
- ``ionisation_mode``, -1 or 1
- ``max_broken_bonds``, Bond breaks
- ``max_water_losses``, Additional water losses
- ``mz_precision``, Relative precision (ppm)

Example response:

.. code-block:: javascript

   {
      "success": true,
      "id": "844bcea5-058b-4b7f-8d29-ba2cc131a568"
   }

Making job public
=================

By default jobs can only be seen by the user that submitted it.
An additional command is needed to make it visible for anyone.

Query file (query.json):

.. code-block:: javascript

   {
      "description": "New description",
      "ms_filename": "New file name for MS data",
      "is_public": true
   }

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar -d @query.json -X PUT http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568

http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568 can now be shared and shown in a web-browser.
When job is not yet completed it will show a status page, after completion the results will be shown.

Poll status
===========

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar http://www.emetabolomics.org/magma/status/844bcea5-058b-4b7f-8d29-ba2cc131a568.json

Where ``844bcea5-058b-4b7f-8d29-ba2cc131a568`` is the job identifier returned by the job submission.

Retry until job has status STOPPED.

Example response:

.. code-block:: javascript

   {
      "status" : "STOPPED",
      "jobid" : "844bcea5-058b-4b7f-8d29-ba2cc131a568"
   }

Fetching results
================

Molecules
---------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/metabolites.json?start=0;limit=10'

Parameters:

- ``start``, Offset in list of molecules
- ``limit``, Maximum nr of molecules to return
- ``scanid``, only return molecules that have hits in scan with this identifier (optional)

Example response:

.. code-block:: javascript

   {
      "totalUnfiltered": 1,
      "total": 1,
      "rows": [{
         "origin": "CHLOROGENIC ACID (1794427)",
         "smiles": "CWVRJTMFETXNAD",
         "probability": 10234.0,
         "molformula": "C16H18O9",
         "assigned": false,
         "reference": "hyperlink ....",
         "mol": "molblock ....",
         "reactionsequence": [],
         "isquery": true,
         "mim": 354.095082,
         "logp": -0.4,
         "level": 1,
         "metid": 1,
         "nhits": 1
      }],
      "scans": [{"rt": null, "id": 1}]}

Long values have been replaced with `....`.

Fields:

- ``metid`` is the molecule identifier.
- ``origin`` is the name of the molecule.

Molecules for one scan
----------------------

To fetch a ranked list of molecules which are annotated for a certain scan.

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/metabolites.json?start=0;limit=10;scanid=1;sort=%5B%7B%22property%22%3A%22score%22%2C%22direction%22%3A%22ASC%22%7D%5D'

``%5B%7B%22property%22%3A%22score%22%2C%22direction%22%3A%22ASC%22%7D%5D``
is the URL encoded (see http://www.faqs.org/rfcs/rfc3986) version of
``[{"property":"score","direction":"ASC"}]`` and orders the molecules with the highest Candidate score first.

Same response as above, but with additional ``score`` and ``deltappm`` fields.

Fragments
---------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/fragments/123/456.json?node=root'

Where ``123`` is the scan identifier and ``456`` is the molecule identifier.

Parameters:

- ``node``, The fragment identifier to fetch children fragments for.

Example response:

.. code-block:: javascript

   {
      "expanded" : true,
      "children" : [
         {
            "deltah" : -1,
            "deltappm" : -0.8824098991817264,
            "mol" : "molblock ....",
            "formula": "C16H17O9",
            "metid" : 23,
            "fragid" : 5,
            "score" : 3,
            "mass" : 370.1263823051,
            "scanid" : 1789,
            "expanded" : true,
            "mz" : 369.119262695312,
            "mslevel" : 1,
            "atoms" : "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15",
            "isAssigned" : false,
            "leaf" : false,
            "children" : [
               {
                  "deltah" : -2,
                  "deltappm" : -1.861685339415437,
                  "mol" : "molblock ....",
                  "formula" : "C7H11O6",
                  "metid" : 23,
                  "fragid" : 6,
                  "score" : 2,
                  "mass" : 115.039519091,
                  "scanid" : 1790,
                  "expanded" : true,
                  "mz" : 113.024360656738,
                  "mslevel" : 2,
                  "atoms" : "14,15,16,20,22,23,24,25",
                  "leaf" : true
               }
            ]
         }
      ]
   }

Fields:

- ``fragid`` is the fragment identifier.
- ``metid`` is the molecule identifier.
- ``scanid`` is the scan identifier.

Chromatogram
------------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/chromatogram.json'

Example response:

.. code-block:: javascript

   {
      "cutoff": 0.0,
      "scans": [{
         "rt": null,
         "ap": 0,
         "intensity": 69989984.0,
         "id": 1
      }]
   }

Fields:

- ``rt`` is the retention time.
- ``ap`` whether scan has molecules assigned to peaks
- ``id`` is the scan identifier.

Mass spectra
------------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/mspectra/1234.json'

Where ``1234`` is the scan identifier.

Example response:

.. code-block:: javascript

   {
      "precursor": {
         "id": 0,
         "mz": 0.0
      },
      "cutoff": 0.0,
      "peaks": [{
         "intensity": 69989984.0,
         "assigned_metid": null,
         "mz": 353.087494
      }],
      "mslevel": 1
   }

Fields:

- ``precursor``, The precursor scan identifier and mz of current scan.
- ``peaks``, list of peaks for current scan.

Extracted ion chromatogram
--------------------------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568/extractedionchromatogram/1234.json'

Where ``1234`` is the molecule identifier.

Example response:

.. code-block:: javascript

   {
      "chromatogram": [{
         "rt": null,
         "intensity": 69989984.0
      }],
      "scans": [{
         "rt": null,
         "id": 1
      }]
   }

Fields:

- ``chromatogram`` is list of intensities of selected molecule for each retention time.
- ``scans`` is list of scan identifiers where selected molecule was fragmented.

Delete job
==========

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar -X DELETE 'http://www.emetabolomics.org/magma/results/844bcea5-058b-4b7f-8d29-ba2cc131a568'
