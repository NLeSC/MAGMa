from pyramid.config import Configurator
from pyramid.events import subscriber, NewRequest

@subscriber(NewRequest)
def extjsurl(event):
    """ Adds extjsroot url to request using extjsroot setting"""
    event.request.extjsroot = event.request.static_url('magmaweb:static/'+event.request.registry.settings['extjsroot'])

def main(global_config, **settings):
    """ This function returns the Magma WSGI application.
    """
    config = Configurator(settings=settings)
    config.add_static_view('static', 'magmaweb:static', cache_max_age=3600)
    config.add_route('home','/')
    config.add_route('results','/results/{jobid}')
    config.add_route('status','/status/{jobid}')
    config.add_route('metabolites.json', '/results/{jobid}/metabolites.json')
    config.add_route('metabolites.csv', '/results/{jobid}/metabolites.csv')
    config.add_route('fragments.json', '/results/{jobid}/fragments/{scanid}/{metid}.json')
    config.add_route('chromatogram.json', '/results/{jobid}/chromatogram.json')
    config.add_route('mspectra.json', '/results/{jobid}/mspectra/{scanid}.json')
    config.add_route('extractedionchromatogram.json','/results/{jobid}/extractedionchromatogram/{metid}.json')
    config.add_route('stderr.txt', '/results/{jobid}/stderr.txt')
    config.add_route('uploaddb', '/uploaddb')
    config.scan('magmaweb')
    return config.make_wsgi_app()
