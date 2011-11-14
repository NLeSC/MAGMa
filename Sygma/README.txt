Sygma README

Requirements
------------

- Extjs installed in Sygma/sygma/static/
- Add config based on production.ini-dist

Run server
----------

To start webserver:

. ../../env/bin/activate
paster serve development.ini --reload
deactivate


To work with models on Python shell
python
from sygma.models import *
from sqlalchemy import create_engine, and_
from sqlalchemy.sql import exists, func
initialize_sql(create_engine('sqlite:///tea_metabolites2_scans_fragments.db'))

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
