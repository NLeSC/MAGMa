from pyramid.config import Configurator
from sqlalchemy import engine_from_config

from sygma.models import initialize_sql

def main(global_config, **settings):
    """ This function returns a Pyramid WSGI application.
    """
    engine = engine_from_config(settings, 'sqlalchemy.')
    initialize_sql(engine)
    config = Configurator(settings=settings)
    config.add_static_view('static', 'sygma:static', cache_max_age=3600)
    config.add_route('home','/')
    config.add_route('metabolites.json', '/metabolites.json')
    config.add_route('fragments.json', '/fragments/{scanid}/{metid}.json')
    config.add_route('chromatogram.json', '/chromatogram.json')
    config.add_route('mspectra.json', '/mspectra/{id}.json')
    config.add_route('scantree.json', '/scan.tree.json')
    config.add_route('chromatogram/hits.json','/chromatogram/hits.json')
    config.add_route('metabolite/scans.json','/metabolite/{id}/scans.json')
    config.scan('sygma')
    return config.make_wsgi_app()

