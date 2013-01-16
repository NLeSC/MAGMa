MAGMaWeb

By Stefan Verhoeven, Lars Ridder and Marijn Sanders.

MAGMaWeb is the web application to start new MAGMa calculations and view the results.

It's a web application with the server-side written in `Python <http://www.python.org>`_ using the `Pyramid web framework <http://www.pylonsproject.org>`_
and the client-side is written in `ExtJS <http://www.sencha.com/products/extjs>`_.

Development installation
========================

1. Download ExtJs at http://sencha.com/products/extjs and unzip in `web/magmaweb/static`.
2. Create users and register jobs see :ref:`User management <user>` or `docs/user.rst <docs/user.rst>`_.
3. Install MAGMa web and it's dependencies

.. code-block:: bash

   virtualenv env
   . env/bin/activate
   cd web
   python setup.py develop

4. Configure and start web server

.. code-block:: bash

   cp production.ini-dist development.ini # make copy of example configuration
   pserve development.ini --reload # starts web server

Goto http://localhost:6543/magma

Production installation
=======================

Additional to the `Development installation`_ to make application more complete/robust/faster.

* `Minimize javscript`_
* Configure reverse http proxy webserver like `nginx`_ to host static content
* Use a faster wsgi python server like `gunicorn`_ or `uWSGI`_
* Install :ref:`job launcher <launcher>` or `docs/launcher.rst <docs/launcher.rst>`_ to perform calculations
* (Optionally) Make ExtJS installation smaller by removing it's `docs`, `builds` directories

nginx
-----

`nginx web server <http://www.nginx.org>`_ is an open source Web server and a reverse proxy server for HTTP,
 SMTP, POP3 and IMAP protocols, with a strong focus on high concurrency, performance and low memory usage.

.. code-block:: bash

   sudo apt-get install nginx-full

Edit /etc/nginx/sites-enabled/default to:

.. code-block:: nginx

   server {
       #listen   80; ## listen for ipv4; this line is default and implied
       #listen   [::]:80 default ipv6only=on; ## listen for ipv6

       server_name $hostname;

       location /magma {
           proxy_set_header        Host $host;
               proxy_set_header        X-Real-IP $remote_addr;
               proxy_set_header        X-Forwarded-For $proxy_add_x_forwarded_for;
               proxy_set_header        X-Forwarded-Proto $scheme;

               client_max_body_size    1000m;
               client_body_buffer_size 128k;
               proxy_connect_timeout   60s;
               proxy_send_timeout      90s;
               proxy_read_timeout      90s;
               proxy_buffering         off;
               proxy_temp_file_write_size 64k;
               proxy_pass http://127.0.0.1:6543;
               proxy_redirect          off;
       }

       location /magma/static/ {
           alias       /home/stefanv/workspace/MAGMa/web/magmaweb/static/;
           expires     30d;
           add_header  Cache-Control public;
           access_log  off;
       }
   }

gunicorn
--------

`Gunicorn wsgi server <http://gunicorn.org/>`_ is a Python WSGI HTTP Server for UNIX.

Edit `development.ini` file by commenting out the `server:main` section with `waitress`.
And remove comment in-front of the `server:main` section with `gunicorn`.

Then start gunicorn with:

.. code-block:: bash

   pip install gunicorn
   pserve development.ini

uWSGI
-----

`uWSGI wsgi server<http://projects.unbit.it/uwsgi/>`_  is a fast,
self-healing and developer/sysadmin-friendly application container server coded in pure C.

The `HttpUwsgiModule <http://wiki.nginx.org/HttpUwsgiModule>`_ is required.

In `production.ini-dist` there is a section for uwsgi configuration.

Change /magma section in /etc/nginx/sites-enabled/default to:

.. code-block:: nginx

    location /magma {
        proxy_set_header        Host $host;
        proxy_set_header        X-Real-IP $remote_addr;
        proxy_set_header        X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header        X-Forwarded-Proto $scheme;

        client_max_body_size    1000m;
        client_body_buffer_size 128k;
        include uwsgi_params;
        uwsgi_pass unix:/tmp/magma.uwsgi.sock;
        uwsgi_param SCRIPT_NAME /magma;
        uwsgi_modifier1 30;
        uwsgi_param  UWSGI_SCHEME   $scheme;
    }

Then start uWSGI with:

.. code-block:: bash

   pip install uwsgi
   uwsgi -H env --ini-paste-logged development.ini

Minimize javascript
-------------------

Install Sencha SDK tools by following instructions at http://www.sencha.com/products/sencha-cmd .

Then concatenate and compress javascript with:

.. code-block:: bash

   cd magmaweb
   sencha build -d static/app -p magmaweb.results-4.1.1a.jsb3
   ln -s magmaweb/static/app/resultsApp-all-4.1.1a.js magmaweb/static/app/resultsApp-all.js

Now not hundreds of seperate javascript files are loaded, but a single javascript file.

Create magmaweb.results.jsb3 file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This only needs to be done if magmaweb.results*.jsb3 does not yet create.

The `sencha create` command does not work for our pages. So we role our own jsb3 writer.

1. Load result page.
2. Goto developers/firebug console
3. Enter `copy(Ext.Loader.history)`
4. Open file `myhistory` and paste clipboard (CTRL-p)
5. Run `perl loader2jsb3.pl myhistory > magmaweb.results-4.1.1a.jsb3`

loader2jsb3.pl looks like:

.. code-block:: perl

   #!/usr/bin/env perl

   use strict;
   use warnings;
   use JSON;

   my %paths = (
      'Ext' => 'static/extjs-4.2.0/src',
      'Ux'  => 'static/extjs-4.2.0/examples/ux',
      'Esc' => 'static/esc',
      'App' => 'static/app'
   );
   my @files;

   while (<>) {
     my $line = $_;
     chomp($line);
     for my $dep (split(/,/,$line)) {
       my ($path, $name) = $dep =~ /(.*)\.(.*)/;
       $name .= '.js';
       $path =~ s/\./\//g;
       $path .= '/';
       if ($path=~/^Esc\/magmaweb/) {
           $path =~ s/^Esc\/magmaweb/$paths{App}/;
       } elsif ($path=~/^Esc/) {
           $path =~ s/^Esc/$paths{Esc}/;
       } elsif ($path=~/^Ext\/ux/) {
           $path =~ s/^Ext\/ux/$paths{Ux}/;
       } else {
   	$path =~ s/^Ext/$paths{Ext}/;
       }
       push(@files, {'path'=> $path, 'name'=> $name});
     }
   }

   print to_json({
     'projectName'=> 'MAGMA web results',
     licenseText=> "Copyright(c) 2011 Netherlands eScience Center",
       "builds"=> [
           {
               "name"=> "All Classes",
               "target"=> "resultsApp-all-4.1.1a.js",
               "compress"=> JSON::true,
               "files"=> \@files
   }
       ],
       "resources"=> []
   }, {pretty=>1});

Running tests
=============

Python
------

Python tests can be run with:

.. code-block:: bash

   pip install nose coverage
   nosetests

Javascript
----------

The ExtJS tests can be run by hosting the `web` directory in a web server (like apache or nginx) and
opening `magmaweb/tests/js/SpecRunner.html` and `magmaweb/tests/js/app/run-tests.html` in a web-browser.

Generate documentation
======================

Python
------

Generate Python documentation with

.. code-block:: bash

   pip install sphinx
   cd docs
   make html
   firefox _build/html/index.html

Javascript
----------

Javascript documentation generation with JSDuck.
See https://github.com/senchalabs/jsduck

.. code-block:: bash

   jsduck magmaweb/static/extjs-4.2.0/src magmaweb/static/extjs-4.2.0/examples/ux \
   magmaweb/static/d3/d3.v3.js magmaweb/static/esc magmaweb/static/app --builtin-classes \
   --output jsdoc --images magmaweb/static/extjs-4.2.0/docs/images
   firefox jsdoc/index.html
