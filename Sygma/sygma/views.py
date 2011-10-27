from sygma.models import DBSession
from sygma.models import Metabolite, Scan, Peak
from pyramid.response import Response
from pyramid.view import view_config
from sqlalchemy.sql.expression import desc, asc
import simplejson as json

@view_config(route_name='home', renderer='templates/results.pt')
def index(request):
    return { 'cutoff': 200000 }

@view_config(route_name='metabolites.json', renderer='json')
def metabolitesjson(request):
    dbsession = DBSession()
    mets = []
    start = int(request.params['start'])
    limit = int(request.params['limit'])
    q = dbsession.query(Metabolite)
#    if (request.params['filter']):

    if ('filter' in request.params):
        for filter in json.loads(request.params['filter']):
            column = Metabolite.__dict__[filter['field']]
            if (filter['type'] == 'numeric'):
                if (filter['comparison'] == 'eq'):
                    q = q.filter(column==filter['value'])
                if (filter['comparison'] == 'gt'):
                    q = q.filter(column>filter['value'])
                if (filter['comparison'] == 'lt'):
                    q = q.filter(column<filter['value'])
            if (filter['type'] == 'string'):
                q = q.filter(column.contains(filter['value']))
            if (filter['type'] == 'list'):
                q = q.filter(column.in_(filter['value']))
            if (filter['type'] == 'boolean'):
                q = q.filter(column==filter['value'])

    total = q.count()

    if ('sort' in request.params):
        for col in json.loads(request.params['sort']):
            if (col['direction'] == 'DESC'):
                q = q.order_by(desc(col['property']))
            elif (col['direction'] == 'ASC'):
                q = q.order_by(asc(col['property']))
    else:
        q = q.order_by(desc(Metabolite.probability), Metabolite.metid)

    # filter on level, probability, reaction seq, formula en scanid en (scanid,mz)
    for met in q[start:(limit+start)]:
        mets.append({
            'id': met.metid,
            'mol': met.mol,
            'level': met.level,
            'probability': met.probability,
            'reactionsequence': met.reactionsequence,
            'smiles': met.smiles,
            'molformula': met.molformula,
            'isquery': met.isquery,
            'origin': met.origin,
            'nhits': met.nhits
        })

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
            'intensity': scan.basepeakintensity,
            'hashit': False
        })

    return scans

@view_config(route_name='mspectra.json', renderer='json')
def mspectrajson(request):
    dbsession = DBSession()
    peaks = []
    for peak in dbsession.query(Peak).filter_by(scanid=request.matchdict['id']):
        peaks.append({
            'mz': peak.mz,
            'intensity': peak.intensity
        })

    return peaks

@view_config(route_name='scantree.json', renderer='json')
def scantree(request):
    return {}
