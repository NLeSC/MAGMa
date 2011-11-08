from sygma.models import DBSession
from sygma.models import Metabolite, Scan, Peak, Fragment, Run
from pyramid.response import Response
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPNotFound
from sqlalchemy.sql.expression import desc, asc
from sqlalchemy.sql import exists, func
from sqlalchemy import and_
from sqlalchemy.orm.exc import NoResultFound
import simplejson as json

@view_config(route_name='home', renderer='templates/home.mak')
def home(request):
    return dict()

@view_config(route_name='results', renderer='results.mak')
def results(request):
    run = DBSession().query(Run).one()
    maxmslevel = DBSession().query(func.max(Scan.mslevel)).scalar()
    return dict(run=run, maxmslevel=maxmslevel)

def extjsgridfilter(q,column,filter):
    if (filter['type'] == 'numeric'):
        if (filter['comparison'] == 'eq'):
            return q.filter(column==filter['value'])
        if (filter['comparison'] == 'gt'):
            return q.filter(column>filter['value'])
        if (filter['comparison'] == 'lt'):
            return q.filter(column<filter['value'])
    elif (filter['type'] == 'string'):
        return q.filter(column.contains(filter['value']))
    elif (filter['type'] == 'list'):
        return q.filter(column.in_(filter['value']))
    elif (filter['type'] == 'boolean'):
        return q.filter(column==filter['value'])

@view_config(route_name='metabolites.json', renderer='json')
def metabolitesjson(request):
    dbsession = DBSession()
    mets = []
    start = int(request.params['start'])
    limit = int(request.params['limit'])
    q = dbsession.query(Metabolite)

    # custom filters
    if ('scanid' in request.params):
        # TODO add score column + order by score
        q = q.join(Fragment.metabolite).filter(Fragment.parentfragid==0).filter(
            Fragment.scanid==request.params['scanid']
        )
    if ('metid' in request.params):
        q = q.filter(Metabolite.metid==request.params['metid'])

    # add nr_scans column
    stmt = DBSession.query(Fragment.metid,func.count('*').label('nr_scans')).filter(
        Fragment.parentfragid==0).group_by(Fragment.metid).subquery()
    q = q.add_column(stmt.c.nr_scans).outerjoin(stmt, Metabolite.metid==stmt.c.metid)

    if ('filter' in request.params):
        for filter in json.loads(request.params['filter']):
            # generic filters
            if (filter['field'] == 'nr_scans'):
                col = stmt.c.nr_scans
            else:
                col = Metabolite.__dict__[filter['field']]
            q = extjsgridfilter(q, col, filter)

    total = q.count()

    if ('sort' in request.params):
        for col in json.loads(request.params['sort']):
            if (col['property'] == 'nr_scans'):
                col2 = stmt.c.nr_scans
            else:
                col2 = Metabolite.__dict__[col['property']]
            if (col['direction'] == 'DESC'):
                q = q.order_by(desc(col2))
            elif (col['direction'] == 'ASC'):
                q = q.order_by(asc(col2))
    else:
        # default sort
        q = q.order_by(desc(Metabolite.probability), Metabolite.metid)

    for (met, nr_scans) in q[start:(limit+start)]:
        r = {
            'metid': met.metid,
            'mol': met.mol,
            'level': met.level,
            'probability': met.probability,
            'reactionsequence': met.reactionsequence,
            'smiles': met.smiles,
            'molformula': met.molformula,
            'isquery': met.isquery,
            'origin': met.origin,
            'nhits': met.nhits,
            'nr_scans': nr_scans
        }
        mets.append(r)

    return { 'total': total, 'rows': mets }

@view_config(route_name='chromatogram.json', renderer='json')
def chromatogramjson(request):
    dbsession = DBSession()
    scans = []
    # TODO add left join to find if scan has metabolite hit
    for scan in dbsession.query(Scan).filter_by(mslevel=1):
        scans.append({
            'id': scan.scanid,
            'rt': scan.rt,
            'intensity': scan.basepeakintensity
        })

    return scans

@view_config(route_name='chromatogram/hits.json', renderer='json')
def chromatogram_hits(request):
    dbsession = DBSession()
    # only metabolite fragment hits
    fq = DBSession.query(Fragment.scanid).filter(Fragment.parentfragid==0)

    if ('metid' in request.params):
        fq.filter(Fragment.metid==request.params['metid'])

    hits = []
    # only scans which are level 1 and have a metabolite hit
    for hit in dbsession.query(Scan.rt,Scan.scanid).filter_by(mslevel=1).filter(Scan.scanid.in_(fq)):
        hits.append({
            'id': hit.scanid,
            'rt': hit.rt,
        })

    return hits

@view_config(route_name='mspectra.json', renderer='json')
def mspectrajson(request):
    dbsession = DBSession()
    scanid = request.matchdict['id']
    peaks = []
    for peak in dbsession.query(Peak).filter_by(scanid=scanid):
        peaks.append({
            'mz': peak.mz,
            'intensity': peak.intensity
        })

    cutoff = DBSession().query(Scan.basepeakintensity*Run.msms_intensity_cutoff).filter(Scan.scanid==scanid).scalar()
    return { 'peaks': peaks, 'cutoff': cutoff }

@view_config(route_name='scantree.json', renderer='json')
def scantree(request):
    return {}

@view_config(route_name='metabolite/scans.json', renderer='json')
def metabolitescans(request):
    metid = request.matchdict['id']
    scans = []
    for frag in DBSession().query(Fragment.scanid).filter(Fragment.metid==metid).filter(Fragment.parentfragid==0):
        scans.append(frag.scanid)
    chromatogram = []
    # fetch avg mz of metabolite fragment
    mzq = DBSession().query(func.avg(Fragment.mz)).filter(Fragment.metid==metid).filter(Fragment.parentfragid==0).scalar()
    if (mzq):
        mzoffset = DBSession().query(Run.mz_precision).scalar()
        # fetch max intensity of peaks with mz = mzq+-mzoffset
        for (rt,intens) in DBSession().query(Scan.rt,func.max(Peak.intensity)).outerjoin(Peak, and_(Peak.scanid==Scan.scanid,Peak.mz.between(mzq-mzoffset,mzq+mzoffset))).filter(Scan.mslevel==1).group_by(Scan.rt).order_by(asc(Scan.rt)):
            chromatogram.append({
                'rt': rt,
                'intensity': intens or 0
            })
    return {'scans': scans, 'chromatogram': chromatogram}

@view_config(route_name='fragments.json', renderer='json')
def fragments(request):
    node = request.params['node']
    def q():
        return DBSession().query(Fragment,Metabolite.mol,Scan.mslevel).join(Metabolite).join(Scan)

    def fragment2json(row):
        (frag, mol, mslevel) = row
        f = {
            'fragid': frag.fragid,
            'scanid': frag.scanid,
            'metid': frag.metid,
            'score': frag.score,
            'mol': mol,
            'atoms': frag.atoms,
            'mz': frag.mz,
            'mass': frag.mass,
            'deltah': frag.deltah,
            'mslevel': mslevel,
        }
        if (len(frag.children)>0):
            f['expanded'] = False
            f['leaf'] = False
        else:
            f['expanded'] = True
            f['leaf'] = True
        return f

    # parent metabolite
    if (node == ''):
        try:
            row = q().filter(
                Fragment.scanid==request.matchdict['scanid']).filter(
                Fragment.metid==request.matchdict['metid']).filter(
                Fragment.parentfragid==0).one()
        except NoResultFound:
            return HTTPNotFound();
        metabolite = fragment2json(row)
        metabolite['children'] = []
        for frow in q().filter(Fragment.parentfragid==metabolite['fragid']):
            metabolite['expanded'] = True
            metabolite['children'].append(fragment2json(frow))
        return { 'children': metabolite, 'expanded': True}
    # fragments
    else:
        fragments = []
        for row in q().filter(Fragment.parentfragid==node):
            fragments.append(fragment2json(row))
        return fragments