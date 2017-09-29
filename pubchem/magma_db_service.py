from waitress import serve
from pyramid.config import Configurator

from magma import PubChemEngine,HmdbEngine


def molecules(request):
    query = request.json_body
    select, incl_halo = query
    print select
    selected_db = request.matchdict['dbname']
    if selected_db == 'hmdb':
        db_engine = HmdbEngine(online = False)
    elif selected_db == 'kegg':
        db_engine = PubChemEngine('kegg', incl_halo=True, online = False)
    elif selected_db == 'pubchem':
        db_engine = PubChemEngine('pubchem', incl_halo=True, online = False)
    return db_engine.query_local(select, incl_halo)

if __name__ == '__main__':

    config = Configurator()
    config.add_route('molecules', '/magma/molecules/{dbname}')
    config.add_view(molecules, route_name='molecules', renderer='json')
    app = config.make_wsgi_app()
    serve(app, host='0.0.0.0', port=8080)

