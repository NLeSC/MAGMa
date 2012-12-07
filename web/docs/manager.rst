.. _manager:

Job manager
===========

To perform MAGMa calculations the web application offloads the job to a manager deamon.

The job manager is in the `jobmanager/` directory of the MAGMa repository.

The example configuration expects a jobmanager deamon running on `http://localhost:9998/job`.

The status of a job managed by the job manager is pushed back to the web application using a callback url.
The jobmanager is authenticated against the web application using it's localhost ip adress.