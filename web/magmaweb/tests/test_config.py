'''
Created on May 31, 2013

@author: verhoes
'''
import unittest
from mock import Mock, call, patch
from pyramid import testing
from pyramid.config import Configurator
from magmaweb.config import configure, jsonhtml_renderer_factory
from magmaweb.user import RootFactory, JobIdFactory


class TestConfigure(unittest.TestCase):

    def setUp(self):
        self.config = Mock(Configurator)
        self.settings = {
                         'jobfactory.root_dir': '/tmp/jobs',
                         'mako.directories': 'magmaweb:templates',
                         'extjsroot': 'ext',
                         'sqlalchemy.url': 'sqlite:///:memory:',
                         'cookie.secret': 'aepeeV6aizaiph5Ae0Reimeequuluwoh',
                         'cookie.path': '/magma',
                         'monitor_user': 'jobmanager',
                         }
        self.config.get_settings.return_value = self.settings

    def testTransaction(self):
        configure(self.config)
        self.config.include.assert_called_once_with('pyramid_tm')

    def testAuth(self):
        configure(self.config)
        self.assertTrue(self.config.set_authentication_policy.called)
        self.assertTrue(self.config.set_authorization_policy.called)
        self.config.set_root_factory.assert_called_once_with(RootFactory)

    def testViews(self):
        configure(self.config)
        c = self.config
        c.add_renderer.assert_called_once_with('jsonhtml',
                                               jsonhtml_renderer_factory)
        c.add_static_view.assert_called_once_with('static',
                                                  'magmaweb:static',
                                                  cache_max_age=3600)
        c.scan.assert_called_once_with('magmaweb', ignore='magmaweb.tests')

    def testRoutes(self):
        configure(self.config)
        calls = [
                 call('home', '/'),
                 call('help', '/help'),
                 call('login', '/login'),
                 call('defaults.json', '/defaults.json'),
                 call('startjob', '/start'),
                 call('jobfromscratch', '/results/'),
                 call('uploaddb', '/uploaddb'),
                 call('workspace', '/workspace'),
                 call('access_token', '/access_token.json'),
                 call('logout', '/logout'),
                 ]
        self.config.add_route.assert_has_calls(calls, any_order=True)

    def testJobRoutes(self):
        configure(self.config)

        def add_job_route(name, pattern):
            return call(name, pattern,
                 traverse='/{jobid}',
                 factory=JobIdFactory)

        calls = [
                 add_job_route('status.json', '/status/{jobid}.json'),
                 add_job_route('status', '/status/{jobid}'),
                 add_job_route('results', '/results/{jobid}'),
                 add_job_route('metabolites.json',
                               '/results/{jobid}/metabolites.json'),
                 add_job_route('metabolites.csv',
                               '/results/{jobid}/metabolites.csv'),
                 add_job_route('metabolites.sdf',
                               '/results/{jobid}/metabolites.sdf'),
                 add_job_route('fragments.json',
                               '/results/{jobid}/fragments/{scanid}/{metid}.json'),
                 add_job_route('chromatogram.json',
                               '/results/{jobid}/chromatogram.json'),
                 add_job_route('mspectra.json',
                               '/results/{jobid}/mspectra/{scanid}.json'),
                 add_job_route('extractedionchromatogram.json',
                               '/results/{jobid}/extractedionchromatogram/{metid}.json'),
                 add_job_route('stderr.txt', '/results/{jobid}/stderr.txt'),
                 add_job_route('runinfo.json',
                               '/results/{jobid}/runinfo.json'),
                 add_job_route('rpc.add_structures',
                               '/rpc/{jobid}/add_structures'),
                 add_job_route('rpc.add_ms_data', '/rpc/{jobid}/add_ms_data'),
                 add_job_route('rpc.metabolize', '/rpc/{jobid}/metabolize'),
                 add_job_route('rpc.metabolize_one',
                               '/rpc/{jobid}/metabolize_one'),
                 add_job_route('rpc.annotate', '/rpc/{jobid}/annotate'),
                 add_job_route('rpc.assign', '/rpc/{jobid}/assign'),
                 add_job_route('rpc.unassign', '/rpc/{jobid}/unassign'),
                 ]
        self.config.add_route.assert_has_calls(calls, any_order=True)

    def testDefaults_None(self):
        configure(self.config)
        calls = [call(auto_register=False), call(restricted=False)]
        self.config.add_settings.assert_has_calls(calls, any_order=True)

    def testDefaults_True(self):
        self.settings['auto_register'] = True
        self.settings['restricted'] = True
        configure(self.config)
        calls = [call(auto_register=True), call(restricted=True)]
        self.config.add_settings.assert_has_calls(calls, any_order=True)

    @patch('magmaweb.config.init_user_db')
    def testDatabase(self, initdb):
        configure(self.config)
        self.assertEquals(str(initdb.call_args[0][0].url),
                          'sqlite:///:memory:')


class Test_jsonhtml_renderer_factory(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()

    def _callFUT(self, name):
        return jsonhtml_renderer_factory(name)

    def test_body(self):
        renderer = self._callFUT(None)
        result = renderer({'a': 1}, {})
        self.assertEqual(result, '{"a": 1}')

    def test_content_type(self):
        request = testing.DummyRequest()
        renderer = self._callFUT(None)
        renderer({'a': 1}, {'request': request})
        self.assertEqual(request.response.content_type, 'text/html')
