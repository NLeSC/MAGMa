jobmanager
==========

Web service to submit jobs via a JavaGAT supported broker.

Submit job to web service using a HTTP POST request.
Reports status of job to submitter using a callback url.

Requirements
------------

- Maven  (http://maven.apache.org)
- JavaGAT (http://www.cs.vu.nl/ibis/)
- Program tarball, to execute on remote system

Install
-------

1. Download/Unzip/Build JavaGAT
1.2 Download from http://gforge.cs.vu.nl/gf/project/javagat/

.. code-block:: bash

   unzip JavaGAT-2.1.3-src.zip
   cd JavaGAT-2.1.3
   ant

2. Add JavaGAT jars in local maven repository, so they can be used as maven dependencies

.. code-block:: bash

   export GAT_VERSION=2.1.3
   cd JavaGAT-$GAT_VERSION
   mvn install:install-file -Dfile=lib/GAT-API.jar -DartifactId=GAT-API -Dversion=$GAT_VERSION -DgroupId=org.gridlab.gat -Dpackaging=jar -DgeneratePom=true
   mvn install:install-file -Dfile=lib/GAT-engine.jar -DartifactId=GAT-engine -Dversion=$GAT_VERSION -DgroupId=org.gridlab.gat -Dpackaging=jar -DgeneratePom=true
   mvn install:install-file -Dfile=lib/GAT-tests.jar -DartifactId=GAT-tests -Dversion=$GAT_VERSION -DgroupId=org.gridlab.gat -Dpackaging=jar -DgeneratePom=true
   mvn install:install-file -Dfile=lib/ibis-util-2.3-pre.jar -DartifactId=ibis-util -Dversion=2.3-pre -DgroupId=ibis -Dpackaging=jar -DgeneratePom=true

3. Build uber-jar or execute from maven.
3.1. Uber-jar, to start on other machine the `magmajobmanager-1.3-SNAPSHOT.jar` file and `adaptors` directory must be copied.

.. code-block:: bash

   mvn package
   java -Dgat.adaptor.path=JavaGAT-2.1.3/lib/adaptors -jar target/magmajobmanager-1.3-SNAPSHOT.jar

3.2 Execute from maven

.. code-block:: bash

   mvn exec:java

During start-up of the web-service you are request to give your password for Glite certificate.
If you enter something it will use Glite to submit/run jobs.
When nothing (just enter) is entered it will use the local machine to run jobs using the JavaGAT localQ broker.

Usage
-----

Web service listens on http://localhost:9998 .

Create a directory with an input file and script

.. code-block:: bash

   mkdir myjob
   cd myjob
   echo 'Lorem ipsum' > input_file
   echo 'hostname;date;wc -l input_file > output_file' > runme.sh

Create a json file (query.json)

.. code-block:: json

   {
      "jobdir": "<absolute path to myjob directory>/",
      "executable": "/bin/sh",
      "prestaged": ["runme.sh", "input_file"],
      "poststaged": ["output_file"],
      "stderr": "stderr.txt",
      "stdout": "stdout.txt",
      "arguments": ["runme.sh"],
      "status_callback_url": "http://localhost/job/myjob/status"
   }

Then submit it

.. code-block:: bash

   curl -H "Content-Type: application/json" -H 'Accept: application/json' -X POST -d @query.json http://localhost:9998/job

After a while `output_file`, `stderr.txt` and `stdout.txt` file appear in `myjob` directory.
"http://localhost/job/myjob/status" will have several PUT HTTP requests send to it.
The PUT requestes contain job statuses like PRE_STAGING, RUNNING, POST_STAGING, STOPPED.

Callback authentication
^^^^^^^^^^^^^^^^^^^^^^^

The status callbacks uses MAC Access Authentication.
The MAC key indentifier and MAC key must be obtained from the provider.

Documentation
-------------

A maven site can be generated with

.. code-block:: bash

   mvn site
   firefox target/site/index.html

