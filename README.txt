MAGMaWeb README

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
----------

1. Minimize js, see chapter below
2. Create logs subdir in MAGMaWeb dir for wsgi server logs
3. Start wsgi server using production.ini
  pserve production.ini
4. Configure reverse proxy webserver like nginx or lighttpd:

Example lighttpd example:

# add expire module before compress module in server.modules

# magmaweb
# for development disable expire
compress.filetype += ( "application/javascript", "applicaton/json" )
alias.url += ( "/magmaweb/static" => "/home/stefanv/workspace/MAGMaWeb/magmaweb/static" )
$HTTP["url"] =~ "/magmaweb" {
  $HTTP["url"] !~ "/magmaweb/static" {
    proxy.server = (
      "/" => (
        "application" => ( "host" => "127.0.0.1", "port" => 6543 )
      )
    )
  }
  $HTTP["url"] =~ "/magmaweb/static" {
    expire.url = ( "" => "access plus 7 days" )
  }
}

Documentation
-------------

Python documentation generation with
  cd docs
  make html

Javascript documentation generation with
  jsduck magmaweb/static/ext-4.0.7/src magmaweb/static/ext-4.0.7/examples/ux magmaweb/static/d3/d3.js magmaweb/static/esc magmaweb/static/app --builtin-classes --output jsdoc --images magmaweb/static/ext-4.0.7/docs/images
or minimal
  jsduck magmaweb/static/esc magmaweb/static/app --output jsdoc


Minimize js
-----------

Create magmaweb.results.jsb3 file
===============================

This only needs to be done if magmaweb.results.jsb3 does not yet create.

export PATH=$PATH:/home/stefanv/SenchaSDKTools-1.2.3/:/home/stefanv/SenchaSDKTools-1.2.3/command/:/home/stefanv/SenchaSDKTools-1.2.3/lib/:home/stefanv/SenchaSDKTools-1.2.3/command/:/home/stefanv/SenchaSDKTools-1.2.3/appbuilder/:/home/stefanv/SenchaSDKTools-1.2.3/jsbuilder/
cd MAGMaWeb/magmaweb
# in results.mak comment resultsApp-all.js, so all dependencies are dynamicly loaded
sencha create jsb -v -p magmaweb.results.jsb3 -a 'http://localhost/magmaweb/results?jobid=2a398f64-3522-11e1-ac1a-0800272c0b38'
# create jsb exits with error, but output file is ok

Edit magmaweb.results.jsb3
- Alter projectName and license
- Remove app-all.js section
- Rename all-classes.js to app/resultsApp-all.js
- Add '"compress": true,' to 'All classes' build

Build magmaweb.results.jsb3
=========================

sencha build -v -d static/app -p magmaweb.results.jsb3
# in results.mak uncomment resultsApp-all.js

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
