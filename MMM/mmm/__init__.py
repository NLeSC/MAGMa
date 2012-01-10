import os.path
from pyramid.config import Configurator
from pyramid_beaker import session_factory_from_settings
from pyramid.events import subscriber, NewRequest

@subscriber(NewRequest)
def extjsurl(event):
    """ Adds extjsroot url to request using extjsroot setting"""
    event.request.extjsroot = event.request.static_url('mmm:static/'+event.request.registry.settings['extjsroot'])

def main(global_config, **settings):
    """ This function returns the MMM WSGI application.
    """
    session_factory = session_factory_from_settings(settings)
    config = Configurator(settings=settings)
    config.set_session_factory(session_factory)
    config.add_static_view('static', 'mmm:static', cache_max_age=3600)
    config.add_route('home','/')
    config.add_route('results','/results')
    config.add_route('metabolites.json', '/metabolites.json')
    config.add_route('fragments.json', '/fragments/{scanid}/{metid}.json')
    config.add_route('chromatogram.json', '/chromatogram.json')
    config.add_route('mspectra.json', '/mspectra/{scanid}.json')
    config.add_route('extractedionchromatogram.json','/extractedionchromatogram/{metid}.json')
    config.scan('mmm')
    return config.make_wsgi_app()
