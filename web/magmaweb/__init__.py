import json
from pyramid.config import Configurator
from pyramid.authorization import ACLAuthorizationPolicy
from pyramid.authentication import AuthTktAuthenticationPolicy
from pyramid_ipauth import IPAuthenticationPolicy
from pyramid_multiauth import MultiAuthenticationPolicy
from sqlalchemy import engine_from_config
from magmaweb.user import init_user_db, RootFactory, JobIdFactory


def main(global_config, **settings):
    """This function returns the Magma WSGI application.
    """
    config = Configurator(settings=settings)

    config.include('pyramid_tm')

    authn_policy1 = AuthTktAuthenticationPolicy(
        secret=settings['cookie.secret'],
        path=settings['cookie.path']
    )
    # TODO Make jobmanager_ipadr a setting
    jobmanager_ipadr = "127.0.*.*"
    authn_policy2 = IPAuthenticationPolicy(jobmanager_ipadr, "jobmanager",
                                           ["monitor"], '127.0.0.1')
    authn_policy = MultiAuthenticationPolicy([authn_policy1, authn_policy2])
    config.set_authentication_policy(authn_policy)
    config.set_authorization_policy(ACLAuthorizationPolicy())
    config.set_root_factory(RootFactory)

    config.add_renderer('jsonhtml', jsonhtml_renderer_factory)

    config.add_static_view('static', 'magmaweb:static', cache_max_age=3600)

    config.add_route('home', '/')  # allow everyone
    config.add_route('defaults.json', '/defaults.json')  # allow everyone

    config.add_route('startjob', '/start')
    config.add_route('jobfromscratch', '/results/')  # calc
    config.add_route('uploaddb', '/uploaddb')  # calc

    config.add_route('status.json', '/status/{jobid}.json')  # my job
    config.add_route('status', '/status/{jobid}')  # my job

    # JobFactory + traverse
    def add_job_route(name, pattern):
        config.add_route(name, pattern, traverse='/{jobid}',
                         factory=JobIdFactory)

    add_job_route('results', '/results/{jobid}')  # my job
    add_job_route('metabolites.json', '/results/{jobid}/metabolites.json')  # my job
    add_job_route('metabolites.csv', '/results/{jobid}/metabolites.csv')  # my job
    add_job_route('metabolites.sdf', '/results/{jobid}/metabolites.sdf')  # my job
    add_job_route('fragments.json',
                  '/results/{jobid}/fragments/{scanid}/{metid}.json')  # my job
    add_job_route('chromatogram.json', '/results/{jobid}/chromatogram.json')  # my job
    add_job_route('mspectra.json', '/results/{jobid}/mspectra/{scanid}.json')  # my job
    add_job_route('extractedionchromatogram.json',
                  '/results/{jobid}/extractedionchromatogram/{metid}.json')  # my job
    add_job_route('stderr.txt', '/results/{jobid}/stderr.txt')  # my job
    add_job_route('runinfo.json', '/results/{jobid}/runinfo.json')  # my job

    add_job_route('rpc.add_structures', '/rpc/{jobid}/add_structures')  # my job + calc
    add_job_route('rpc.add_ms_data', '/rpc/{jobid}/add_ms_data')  # my job + calc
    add_job_route('rpc.metabolize', '/rpc/{jobid}/metabolize')  # my job + calc
    add_job_route('rpc.metabolize_one', '/rpc/{jobid}/metabolize_one')  # my job + calc
    add_job_route('rpc.annotate', '/rpc/{jobid}/annotate')  # my job + calc
    add_job_route('rpc.set_description', '/rpc/{jobid}/set_description')  # my job + calc
    add_job_route('rpc.assign', '/rpc/{jobid}/assign')  # my job + assigner
    add_job_route('rpc.unassign', '/rpc/{jobid}/unassign')  # my job + assigner

    # user preferences (email, displayName) and
    # list of jobs and their history
    config.add_route('workspace', '/workspace')  # my jobs logged in

    config.add_route('login', '/login')
    config.add_route('logout', '/logout')

    # Setup connection to user database
    engine = engine_from_config(settings)
    init_user_db(engine)

    config.scan('magmaweb')
    return config.make_wsgi_app()


def jsonhtml_renderer_factory(info):
    """Json renderer with text/html content type

    ExtJS form with file upload requires
    json response with text/html content type.
    See http://extjs.com/deploy/dev/docs/?class=Ext.form.BasicForm hasUpload().
    """
    def _render(value, system):
        request = system.get('request')
        if request is not None:
            response = request.response
            response.content_type = 'text/html'
        return json.dumps(value)
    return _render
