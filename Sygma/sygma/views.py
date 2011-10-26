from sygma.models import DBSession
from sygma.models import Metabolite
from pyramid.response import Response
from pyramid.view import view_config

@view_config(route_name='home', renderer='templates/results.pt')
def index(request):
    return { 'cutoff': 200000 }

@view_config(route_name='metabolites.json', renderer='json')
def metabolitesjson(request):
    dbsession = DBSession()
    mets = []
    for met in dbsession.query(Metabolite):
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

    return { 'total': len(mets), 'rows': mets }
