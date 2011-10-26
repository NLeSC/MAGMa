from sygma.models import DBSession
from sygma.models import Metabolite, Scan, Peak
from pyramid.response import Response
from pyramid.view import view_config

@view_config(route_name='home', renderer='templates/results.pt')
def index(request):
    return { 'cutoff': 200000 }

@view_config(route_name='metabolites.json', renderer='json')
def metabolitesjson(request):
    dbsession = DBSession()
    mets = []
    start = int(request.params['start'])
    limit = int(request.params['limit'])
    for met in dbsession.query(Metabolite)[start:(limit+start)]:
        mets.append({
            'id': met.id,
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

    total = dbsession.query(Metabolite).count()
    return { 'total': total, 'rows': mets }

@view_config(route_name='chromatogram.json', renderer='json')
def chromatogramjson(request):
    dbsession = DBSession()
    scans = []
    # TODO add left join to find if scan has metabolite hit
    for scan in dbsession.query(Scan).filter_by(level=1):
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

