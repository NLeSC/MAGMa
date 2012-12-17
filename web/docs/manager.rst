.. _manager:

Job manager
===========

To perform MAGMa calculations the web application offloads the job to a manager deamon.

The job manager is in the `jobmanager/` directory of the MAGMa repository.

The example configuration expects a jobmanager deamon running on `http://localhost:9998/job`.

The status of a job managed by the job manager is pushed back to the web application using a callback url.

The jobmanager is authenticated against the web application using MAC Access Authentication.
See http://tools.ietf.org/html/draft-ietf-oauth-v2-http-mac-01 .

Installation
------------

1. Create a user for the jobmanager, see User management <user>.
2. Set username as 'monitor_user' key in web config file, so that user is authorized with status update rights.
3. On workspace page generate an access token.
4. Configure jobmanager with supplied token.
