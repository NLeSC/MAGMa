===============
MAGMaWeb README
===============

Requirements
------------

- Extjs installed in MAGMaWeb/magmaweb/static/
- Add config based on production.ini-dist

Run server
----------

To start webserver:

. ../../env/bin/activate
paster serve development.ini --reload
deactivate


To work with models on Python shell
python
from magmaweb.models import *
from sqlalchemy import create_engine, and_
from sqlalchemy.sql import exists, func
initialize_sql(create_engine('sqlite:///tea_metabolites2_scans_fragments.db'))

Production deployment
---------------------

1. Minimize js, see chapter below
2. Create logs subdir in MAGMaWeb dir for wsgi server logs
3. Start wsgi server using production.ini
  pserve production.ini
4. Configure reverse proxy webserver like nginx or lighttpd:

Example lighttpd example:
~~~~~~~~~~~~~~~~~~~~~~~~~

# add expire module before compress module in server.modules

# magmaweb
# for development disable expire
compress.filetype += ( "application/javascript", "applicaton/json" )
alias.url += ( "/magma/static" => "/home/stefanv/workspace/MAGMaWeb/magmaweb/static" )
$HTTP["url"] =~ "/magma" {
  $HTTP["url"] !~ "/magma/static" {
    proxy.server = (
      "/" => (
        "application" => ( "host" => "127.0.0.1", "port" => 6543 )
      )
    )
  }
  $HTTP["url"] =~ "/magma/static" {
    expire.url = ( "" => "access plus 70 days" )
  }
}

For https add setenv module and add
setenv.add-request-header = (
      "X-FORWARDED-PROTOCOL" => "ssl"
    )
Above proxy.serve = ...

Example nginx (http proxy):
~~~~~~~~~~~~~~~~~~~~~~~~~~~

sudo apt-get install nginx-full

Edit /etc/nginx/sites-enabled/default to:
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
        alias       /home/stefanv/workspace/MAGMaWeb/magmaweb/static/;
        expires     30d;
        add_header  Cache-Control public;
        access_log  off;
    }
}

Example nginx+uwsgi:
~~~~~~~~~~~~~~~~~~~~

Add to ini file:
[uwsgi]
socket = /tmp/magma.uwsgi.sock
virtualenv = /home/stefanv/workspace/MAGMaWeb/env64
master = true
processes = 40
chmod = 666

Change /magma sectio in /etc/nginx/sites-enabled/default to:
    location /magma {
        proxy_set_header        Host $host;
            proxy_set_header        X-Real-IP $remote_addr;
            proxy_set_header        X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header        X-Forwarded-Proto $scheme;

            client_max_body_size    1000m;
            client_body_buffer_size 128k;
            include uwsgi_params;
            uwsgi_pass unix:/tmp/magma.uwsgi.sock;
            uwsgi_param SCRIPT_NAME /;
    }

pip install uwsgi
uwsgi --ini-paste-logged development.ini

To protect /magma behind basic authentication use:
        auth_basic "Magma login";
        auth_basic_user_file magma.users;
        uwsgi_param REMOTE_USER $remote_user;
        uwsgi_param AUTH_TYPE Basic;

Password file can be managed with 'htpasswd'.

Documentation
-------------

Python documentation generation with
  cd docs
  make html

Javascript documentation generation with
  jsduck magmaweb/static/ext-4.1.1a/src magmaweb/static/ext-4.1.1a/examples/ux magmaweb/static/d3/d3.v2.js magmaweb/static/esc magmaweb/static/app --builtin-classes --output jsdoc --images magmaweb/static/ext-4.1.1a/docs/images
or minimal
  jsduck magmaweb/static/esc magmaweb/static/app --output jsdoc


Minimize js
-----------

Install Sencha SDK tools.

Create magmaweb.results.jsb3 file
===============================

This only needs to be done if magmaweb.results.jsb3 does not yet create.

The `sencha create` command does not work for our pages. So we role our own jsb3 writer.

1. Load result page.
2. Goto developers/firebug console
3. Enter `copy(Ext.Loader.history)`
4. Open file `myhistory` and paste clipboard (CTRL-p)
5. Run `perl loader2jsb3.pl myhistory > magmaweb.results-4.1.1a.jsb3`

loader2jsb3.pl looks like:
<pre>
#!/usr/bin/env perl

use strict;
use warnings;
use JSON;

my %paths = (
   'Ext' => 'static/ext-4.1.1a/src',
   'Ux'  => 'static/ext-4.1.1a/examples/ux',
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
</pre>

Build magmaweb.results.jsb3
===========================

sencha build -v -d static/app -p magmaweb.results.jsb3
# in results.mak uncomment resultsApp-all.js

User management
---------------

Users can be added/edited from python prompt.

Start with:
python
from transaction import commit
from sqlalchemy import create_engine
from magmaweb.user import init_user_db, User
engine = create_engine('sqlite:///data/users.db')
init_user_db(engine)

Add user
========

user = User(u'stefanv2', u'Stefan Verhoeven',
            u's.verhoeven@esciencecenter.nl', 'mypassword')
User.add(user)
commit()

Fetch user
==========

user = User.by_id('stefanv2')

Update user
===========

user.displayname = 'Stefan second account'
User.add(user)
commit()

Change password
===============

user.password = 'mypw'
User.add(user)
commit()

SQL cookbook
------------

print DBSession.query(Metabolite).first()

Rt in chromatogram where metid=352 has hits
SELECT scanid, rt FROM scans WHERE mslevel=1 AND scanid IN (SELECT scanid FROM fragments WHERE metid=352 AND parentfragid=0);
DBSession.query(Scan.rt,Scan.scanid).filter_by(mslevel=1).filter(Scan.scanid.in_(DBSession.query(Fragment.scanid).filter(Fragment.parentfragid==0))).all()

Rt in chromatogram where any metabolites have hits
SELECT scanid, rt FROM scans WHERE mslevel=1 AND scanid IN (SELECT scanid FROM fragments WHERE parentfragid=0);
print DBSession.query(Scan.rt,Scan.scanid).filter_by(mslevel=1).filter(Scan.scanid.in_(DBSession.query(Fragment.scanid).filter(Fragment.parentfragid==0))).all()

Select metabolites with score on scanid=2682
SELECT m.*,score FROM metabolites m JOIN fragments USING (metid) WHERE scanid=2682 AND parentfragid=0;
print DBSession.query(Metabolite,Fragment.score).join(Fragment).filter(Fragment.parentfragid==0).filter(Fragment.scanid==2682)

Select metabolites with nr_scans:
SELECT m.*,nr_scans FROM metabolites m
LEFT JOIN (SELECT metid, count(*) nr_scans FROM fragments WHERE parentfragid=0 GROUP BY metid) USING (metid)
;
SELECT metabolites.metid AS metabolites_metid, anon_1.nr_scans AS anon_1_nr_scans
FROM metabolites LEFT OUTER JOIN (SELECT fragments.metid AS metid, count(*) AS nr_scans
FROM fragments
WHERE fragments.parentfragid = 0 GROUP BY fragments.metid) AS anon_1 ON metabolites.metid = anon_1.metid;

from sqlalchemy.sql import func
stmt = DBSession.query(Fragment.metid,func.count('*').label('nr_scans')).filter(Fragment.parentfragid==0).group_by(Fragment.metid).subquery()
print DBSession.query(Metabolite.metid, stmt.c.nr_scans).outerjoin(stmt, Metabolite.metid==stmt.c.metid)

Select fragments on scanid=2682
SELECT * FROM fragments WHERE scanid=2682;
print DBSession.query(Fragment).filter(Fragment.scanid==2682)

Select child scans of scanid=2682
SELECT * FROM scans WHERE precursorscanid=2682;

# fragment with mol and with haschildren column
SELECT f.*,m.mol,
EXISTS(SELECT * FROM fragments WHERE parentfragid=f.fragid) haschildren
FROM fragments f JOIN metabolites m USING (metid) WHERE  scanid=1263 AND metid=352 AND parentfragid=0
;
r = DBSession.query(Fragment).filter(Fragment.scanid==1263).filter(Fragment.metid==352).filter(Fragment.parentfragid==0).first()
len(r.chilren)


# chromatogram of a metabolite based on its avg mz (from fragments table) +- 0.005
SELECT
s.rt,
p.intensity
FROM
scans s
LEFT JOIN peaks p ON p.scanid = s.scanid
AND p.mz BETWEEN
(SELECT avg(mz)-0.005 FROM fragments WHERE metid = 114 AND parentfragid = 0)
AND
(SELECT avg(mz)+0.005 FROM fragments WHERE metid = 114 AND parentfragid = 0)
WHERE
mslevel = 1
;
SELECT rt, max(p.intensity)
FROM scans s
LEFT JOIN peaks p ON p.scanid=s.scanid AND p.mz BETWEEN
207.0613350266667 AND 207.0713350266667
WHERE mslevel=1
GROUP BY rt ORDER BY rt
;
mzq = DBSession().query(func.avg(Fragment.mz)).filter(Fragment.metid==114).filter(Fragment.parentfragid==0).one()[0]
mzoffset = 0.005
print DBSession().query(Scan.rt,func.max(Peak.intensity)).outerjoin(Peak,and_(Peak.scanid==Scan.scanid,Peak.mz.between(mzq-mzoffset,mzq+mzoffset))).group_by(Scan.rt).order_by(Scan.rt)

CREATE UNIQUE INDEX IF NOT EXISTS scanmz ON peaks (scanid,mz);
CREATE UNIQUE INDEX IF NOT EXISTS scanpk ON scans (scanid);
CREATE INDEX IF NOT EXISTS scanlvl ON scans (mslevel);

SELECT basepeakintensity*msms_intensity_cutoff FROM scans,run WHERE scanid=1104;
print DBSession().query(Scan.basepeakintensity*Run.msms_intensity_cutoff).filter(Scan.scanid==1104)

# fetch precursors of a scan
print DBSession().query(Scan.scanid, Scan.mslevel, Scan.precursormz, Scan.precursorscanid).filter(Scan.scanid==64)
t = DBSession().query(Scan).filter(Scan.scanid==64).one()
t.precursor
