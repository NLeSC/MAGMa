import json
from pyramid.config import Configurator
from pyramid.events import subscriber, NewRequest

@subscriber(NewRequest)
def extjsurl(event):
    """Adds extjsroot url to request using extjsroot setting"""
    event.request.extjsroot = event.request.static_url('magmaweb:static/'+event.request.registry.settings['extjsroot'])

def main(global_config, **settings):
    """This function returns the Magma WSGI application.
    """
    config = Configurator(settings=settings)
    config.add_renderer('jsonhtml', jsonhtml_renderer_factory)

    config.add_static_view('static', 'magmaweb:static', cache_max_age=3600)
    config.add_route('home','/')
    config.add_route('jobfromscratch', '/results/')
    config.add_route('results','/results/{jobid}')
    config.add_route('status.json','/status/{jobid}.json')
    config.add_route('status','/status/{jobid}')
    config.add_route('metabolites.json', '/results/{jobid}/metabolites.json')
    config.add_route('metabolites.csv', '/results/{jobid}/metabolites.csv')
    config.add_route('metabolites.sdf', '/results/{jobid}/metabolites.sdf')
    config.add_route('fragments.json', '/results/{jobid}/fragments/{scanid}/{metid}.json')
    config.add_route('chromatogram.json', '/results/{jobid}/chromatogram.json')
    config.add_route('mspectra.json', '/results/{jobid}/mspectra/{scanid}.json')
    config.add_route('extractedionchromatogram.json','/results/{jobid}/extractedionchromatogram/{metid}.json')
    config.add_route('stderr.txt', '/results/{jobid}/stderr.txt')
    config.add_route('uploaddb', '/uploaddb')
    config.add_route('runinfo.json', '/results/{jobid}/runinfo.json')
    config.add_route('defaults.json', '/defaults.json')

    config.add_route('rpc.add_structures', '/rpc/{jobid}/add_structures')
    config.add_route('rpc.add_ms_data', '/rpc/{jobid}/add_ms_data')
    config.add_route('rpc.metabolize', '/rpc/{jobid}/metabolize')
    config.add_route('rpc.metabolize_one', '/rpc/{jobid}/metabolize_one')
    config.add_route('rpc.annotate', '/rpc/{jobid}/annotate')
    config.add_route('rpc.allinone', '/rpc/{jobid}/allinone')
    config.add_route('rpc.set_description', '/rpc/{jobid}/set_description')

    config.scan('magmaweb')
    return config.make_wsgi_app()

def jsonhtml_renderer_factory(info):
    """Json renderer with text/html content type

    ExtJS form with file upload requires json response with text/html content type.

    See http://extjs.com/deploy/dev/docs/?class=Ext.form.BasicForm hasUpload().
    """
    def _render(value, system):
        request = system.get('request')
        if request is not None:
            response = request.response
            response.content_type = 'text/html'
        return json.dumps(value)
    return _render
