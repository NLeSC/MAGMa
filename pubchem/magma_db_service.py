from waitress import serve
from pyramid.config import Configurator

from magma import PubChemEngine,KeggEngine,HmdbEngine


def masses(request):
    db_engines = {
                 'pubchem': PubChemEngine,
                 'kegg': KeggEngine,
                 'hmdb': HmdbEngine
                 }
    queries = request.json_body
    selected_db = request.matchdict['dbname']
    if selected_db == 'hmdb':
        db_engine = HmdbEngine(online = False)
    elif selected_db == 'kegg':
        db_engine = PubChemEngine('kegg', online = False)
    elif selected_db == 'pubchem':
        db_engine = PubChemEngine('pubchem', online = False)
    results = []
    for low, high, charge in queries:
        results += db_engine.query_local(low, high, charge)
    print results
    return results

if __name__ == '__main__':

    config = Configurator()
    config.add_route('masses', '/magma/masses/{dbname}')
    config.add_view(masses, route_name='masses', renderer='json')
    app = config.make_wsgi_app()
    serve(app, host='0.0.0.0', port=8080)