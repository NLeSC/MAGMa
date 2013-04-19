====================
Web service consumer
====================

MAGMa has a nice web user interface.
The user interface recieves/send data using a web service.
The web service can also be used by other applications.
Below is an example how to use the web service using ``curl``, a command line HTTP client.

Submit job
==========

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar -F ms_data_file=@input.mtree \
   -F ms_data_format=mass_tree -F structure_database=pubchem \
   http://www.emetabolomics.org/magma/start

Parameters:

- ``ms_data_file``, file upload of mass spectra data
- ``ms_data_format``, which format ``ms_data_file`` is in. See http://www.emetabolomics.org/magma/help for examples.
- ``structure_database``, Retrieve molecules from a database. Can be pubchem, kegg, hmdb.

TODO parameters incomplete add more to make runnable

Example response:

.. code-block:: javascript

   {
      "success": true,
      "id": "b1eee101-dcc6-435e-baa8-d35e688c408e"
   }

Poll status
===========

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar http://www.emetabolomics.org/magma/status/b1eee101-dcc6-435e-baa8-d35e688c408e.json

Where ``b1eee101-dcc6-435e-baa8-d35e688c408e`` is the job identifier returned by the job submission.

Retry until job has status STOPPED.

Example response:

.. code-block:: javascript

   {
      "status" : "STOPPED",
      "jobid" : "b1eee101-dcc6-435e-baa8-d35e688c408e"
   }

Fetching results
================

Molecules
---------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/b1eee101-dcc6-435e-baa8-d35e688c408e/metabolites.json?start=0;limit=10'

Parameters:

- ``start``, Offset in list of molecules
- ``limit``, Maximum nr of molecules to return
- ``scanid``, only return molecules that have hits in scan with this identifier (optional)

Fragments
---------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/b1eee101-dcc6-435e-baa8-d35e688c408e/fragments.json?scanid=1;metid=2'

Parameters:

- ``scanid``, Fragments on scan with this identifier
- ``metid``, Fragments of metabolite with this identifier
- ``node``, The fragment identifier to fetch children fragments for (optional)

Chromatogram
------------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/b1eee101-dcc6-435e-baa8-d35e688c408e/chromatogram.json'

Mass spectra
------------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/b1eee101-dcc6-435e-baa8-d35e688c408e/mspectra/1234.json'

Where ``1234`` is the scan identifier.

Extracted ion chromatogram
--------------------------

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar 'http://www.emetabolomics.org/magma/results/b1eee101-dcc6-435e-baa8-d35e688c408e/extractedionchromatogram/1234.json'

Where ``1234`` is the molecule identifier.

Delete job
==========

.. code-block:: bash

   curl -c cookie.jar -b cookie.jar -X DELETE 'http://www.emetabolomics.org/magma/results/b1eee101-dcc6-435e-baa8-d35e688c408e'



