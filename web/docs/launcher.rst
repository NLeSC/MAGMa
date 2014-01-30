.. _launcher:

Job launcher
============

To perform MAGMa calculations the web application offloads the job to a launcher deamon.

The job launcher is in the `joblauncher/` directory of the MAGMa repository.

The example configuration expects a joblauncher deamon running on `http://localhost:9998/job`.

The status of a job managed by the job launcher is pushed back to the web application using a callback url.

The joblauncher is authenticated against the web application using MAC Access Authentication.
See http://tools.ietf.org/html/draft-ietf-oauth-v2-http-mac-01 .

Installation
------------

1. Create a user for the joblauncher, see :ref:`User management <user>` or `user.rst <user.rst>`_.
2. Set username as 'monitor_user' key in web config file, so that user is authorized with to update statuses of jobs.
3. Login as the joblauncher user on web site
4. On workspace page generate an access token. Example token object:

    {
        "acesss_token": "bGFpYmFldGhlZXRoZWV0aGFwaGVpZm9vZmFlc2hvcmVpbW9oamluZ2lheG9jaGV6b3VwZW92YWVzaGVhYmFobmdvb3F1ZWlib2thaG5nZWV0ZWVwaG9odGhldXR=",
        "mac_key": "dW9oYml1cGllZ2hhaWNhdWZvaAo=",
        "expires_in": 30758400, # a year
        "token_type": "mac",
        "mac_algorithm": "hmac-sha-1"
    }

5. Configure joblauncher with supplied token.
* Use "access_token" from token object in joblauncher config as "id"
* Use "mac_key" from token object in joblauncher config as "key"
* Use "http(s):\\host:port" where MAGMa web is running at as "scope" in joblauncher config. Examples:
  * if MAGMaWeb is running at http://localhost:6543/magma use http://localhost:6543 as scope
  * if MAGMaWeb is running at https://www.emetabolomics.org/magma use https://www.emetabolomics.org/magma as scope

The joblauncher config will look like:

   macs:
   - id: bGFpYmFldGhlZXRoZWV0aGFwaGVpZm9vZmFlc2hvcmVpbW9oamluZ2lheG9jaGV6b3VwZW92YWVzaGVhYmFobmdvb3F1ZWlib2thaG5nZWV0ZWVwaG9odGhldXR=
     key: dW9oYml1cGllZ2hhaWNhdWZvaAo=
     scope: https://www.emetabolomics.org


