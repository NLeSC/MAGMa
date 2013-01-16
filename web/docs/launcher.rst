.. _launcher:

Job launcher
===========

To perform MAGMa calculations the web application offloads the job to a launcher deamon.

The job launcher is in the `joblauncher/` directory of the MAGMa repository.

The example configuration expects a joblauncher deamon running on `http://localhost:9998/job`.

The status of a job managed by the job launcher is pushed back to the web application using a callback url.

The joblauncher is authenticated against the web application using MAC Access Authentication.
See http://tools.ietf.org/html/draft-ietf-oauth-v2-http-mac-01 .

Installation
------------

1. Create a user for the joblauncher, see :ref:`User management <user>` or `user.rst <user.rst>`_.
2. Set username as 'monitor_user' key in web config file, so that user is authorized with status update rights.
3. On workspace page generate an access token.
4. Configure joblauncher with supplied token.
