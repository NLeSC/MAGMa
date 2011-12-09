from pyramid.config import Configurator
from sqlalchemy import engine_from_config
from pyramid_beaker import session_factory_from_settings

from sygma.models import initialize_sql

def main(global_config, **settings):
    """ This function returns the MSygma WSGI application.
    """
    engine = engine_from_config(settings, 'sqlalchemy.')
    session_factory = session_factory_from_settings(settings)
    initialize_sql(engine)
    config = Configurator(settings=settings)
    config.set_session_factory(session_factory)
    config.add_static_view('static', 'sygma:static', cache_max_age=3600)
    config.add_route('home','/')
    config.add_route('results','/results')
    config.add_route('metabolites.json', '/metabolites.json')
    config.add_route('fragments.json', '/fragments/{scanid}/{metid}.json')
    config.add_route('chromatogram.json', '/chromatogram.json')
    config.add_route('mspectra.json', '/mspectra/{scanid}.json')
    config.add_route('extractedionchromatogram.json','/extractedionchromatogram/{metid}.json')
    config.scan('sygma')
    return config.make_wsgi_app()

