'''
Created on May 31, 2013

@author: verhoes
'''
import json
from pyramid.authorization import ACLAuthorizationPolicy
from pyramid.authentication import AuthTktAuthenticationPolicy
from pyramid.settings import asbool
from pyramid_macauth import MACAuthenticationPolicy
from pyramid_multiauth import MultiAuthenticationPolicy
from sqlalchemy import engine_from_config
from magmaweb.user import init_user_db, RootFactory, JobIdFactory
from magmaweb.user import groupfinder


def configure(config):
    """Configures the MAGMa web app.

    `config` is a instance of :class:Configurator.
    """
    settings = config.get_settings()
    config.include('pyramid_mako')
    config.include('pyramid_tm')

    # for human users
    authn_policy1 = AuthTktAuthenticationPolicy(
        secret=settings['cookie.secret'],
        path=settings['cookie.path'],
        hashalg='sha512',
        callback=groupfinder,
        cookie_name=settings.get('cookie.name', 'auth_tkt'),
        wild_domain=False,
        http_only=True,
    )

    # for service consumers
    # See http://www.rfk.id.au/blog/entry/securing-pyramid-persona-macauth/
    authn_policy2 = MACAuthenticationPolicy.from_settings(settings)
    authn_policy2.find_groups = groupfinder

    auth_policies = [authn_policy1, authn_policy2, ]
    authn_policy = MultiAuthenticationPolicy(auth_policies)
    config.set_authentication_policy(authn_policy)
    config.set_authorization_policy(ACLAuthorizationPolicy())
    config.set_root_factory(RootFactory)

    config.add_renderer('jsonhtml', jsonhtml_renderer_factory)

    config.add_static_view('static', 'magmaweb:static', cache_max_age=3600)

    # for everyone
    config.add_route('home', '/')
    config.add_route('help', '/help')
    config.add_route('login', '/login')

    # for authenticated users
    config.add_route('defaults.json', '/defaults.json')
    config.add_route('startjob', '/start')
    config.add_route('jobfromscratch', '/results/')
    config.add_route('uploaddb', '/uploaddb')
    config.add_route('workspace', '/workspace')
    config.add_route('access_token', '/access_token.json')
    config.add_route('logout', '/logout')

    # JobFactory + traverse
    def add_job_route(name, pattern):
        """"Add route with :class:Job as context"""
        config.add_route(name, pattern, traverse='/{jobid}',
                         factory=JobIdFactory)

    # for job owner
    add_job_route('status.json', '/status/{jobid}.json')
    add_job_route('status', '/status/{jobid}')

    # for authenticated users
    add_job_route('results', '/results/{jobid}')
    add_job_route('molecules.json', '/results/{jobid}/molecules.json')
    add_job_route('molecules.csv', '/results/{jobid}/molecules.csv')
    add_job_route('molecules.sdf', '/results/{jobid}/molecules.sdf')
    add_job_route('fragments.json',
                  '/results/{jobid}/fragments/{scanid}/{molid}.json')
    add_job_route('chromatogram.json', '/results/{jobid}/chromatogram.json')
    add_job_route('mspectra.json', '/results/{jobid}/mspectra/{scanid}.json')
    add_job_route('extractedionchromatogram.json',
                  '/results/{jobid}/extractedionchromatogram/{molid}.json')
    add_job_route('stderr.txt', '/results/{jobid}/stderr.txt')
    add_job_route('stdout.txt', '/results/{jobid}/stdout.txt')
    add_job_route('runinfo.json', '/results/{jobid}/runinfo.json')

    # for job owner
    add_job_route('rpc.add_structures', '/rpc/{jobid}/add_structures')
    add_job_route('rpc.add_ms_data', '/rpc/{jobid}/add_ms_data')
    add_job_route('rpc.metabolize', '/rpc/{jobid}/metabolize')
    add_job_route('rpc.metabolize_one', '/rpc/{jobid}/metabolize_one')
    add_job_route('rpc.annotate', '/rpc/{jobid}/annotate')
    add_job_route('rpc.assign', '/rpc/{jobid}/assign')
    add_job_route('rpc.unassign', '/rpc/{jobid}/unassign')

    # find view_config decorations
    config.scan('magmaweb', ignore='magmaweb.tests')

    # add config defaults and
    # cast config parameter to boolean
    auto_register = asbool(settings.get('auto_register', False))
    config.add_settings(auto_register=auto_register)
    restricted = asbool(settings.get('restricted', False))
    config.add_settings(restricted=restricted)
    ncpus = settings.get('ncpus', 1)
    config.add_settings(ncpus=ncpus)

    # Setup connection to user database
    engine = engine_from_config(settings)
    init_user_db(engine)


def jsonhtml_renderer_factory(info):
    """Json renderer with text/html content type

    ExtJS form with file upload requires
    json response with text/html content type.
    See http://extjs.com/deploy/dev/docs/?class=Ext.form.BasicForm hasUpload().
    """
    def _render(value, system):
        """Return jsonified 'value' with 'text/html' as content type"""
        request = system.get('request')
        if request is not None:
            response = request.response
            response.content_type = 'text/html'
        return json.dumps(value)
    return _render
