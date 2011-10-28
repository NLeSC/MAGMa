Sygma README

To start webserver:

. ../../env/bin/activate
paster serve development.ini --reload
deactivate


To work with models on Python shell
python
from sygma.models import *
from sqlalchemy import create_engine
initialize_sql(create_engine('sqlite:///tea_metabolites2_scans_fragments.db'))
print DBSession.query(Metabolite).first()


SQL cookbook
------------

Rt in chromatogram where metid=352 has hits
SELECT scanid, rt FROM scans WHERE mslevel=1 AND scanid IN (SELECT scanid FROM fragments WHERE metid=352 AND parentfragid=0);
DBSession.query(Scan.rt,Scan.scanid).filter_by(mslevel=1).filter(Scan.scanid.in_(DBSession.query(Fragment.scanid).filter(Fragment.parentfragid==0))).all()

Rt in chromatogram where any metabolites have hits
SELECT scanid, rt FROM scans WHERE mslevel=1 AND scanid IN (SELECT scanid FROM fragments WHERE parentfragid=0);
print DBSession.query(Scan.rt,Scan.scanid).filter_by(mslevel=1).filter(Scan.scanid.in_(DBSession.query(Fragment.scanid).filter(Fragment.parentfragid==0))).all()

Select metabolites with score on scanid=2682
SELECT m.*,score FROM metabolites m JOIN fragments USING (metid) WHERE scanid=2682 AND parentfragid=0;
print DBSession.query(Metabolite,Fragment.score).join(Fragment).filter(Fragment.parentfragid==0).filter(Fragment.scanid==2682)

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