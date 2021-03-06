import unittest
import datetime
import json
from pyramid import testing
from mock import Mock, patch
from magmaweb.views import Views, JobViews, InCompleteJobViews
from magmaweb.job import JobFactory, Job, JobDb, JobQuery
from magmaweb.job import JobError, JobIncomplete
from magmaweb.user import User, JobMeta
try:  # to stay python 2 compatible
    from StringIO import StringIO
except:
    from io import StringIO


class AbstractViewsTestCase(unittest.TestCase):

    def setUp(self):
        self.settings = {'jobfactory.root_dir': '/somedir',
                         'access_token.expires_in': 360,
                         'auto_register': False,
                         'restricted': False,
                         'ncpus': 1}
        self.config = testing.setUp(settings=self.settings)
        self.config.add_route('status.json', '/status/{jobid}.json')

    def tearDown(self):
        testing.tearDown()

    def fake_job(self):
        job = Mock(Job)
        job.id = 'foo'
        job.description = ""
        job.ms_filename = ""
        job.dir = '/somedir/foo'
        job.state = 'STOPPED'
        job.db = Mock(JobDb)
        job.db.runInfo.return_value = 'bla'
        job.db.maxMSLevel.return_value = 3
        job.db.molecules.return_value = {'total': 3, 'rows': [1, 2, 3]}
        job.db.moleculesTotalCount.return_value = 3
        job.db.scansWithMolecules.return_value = [4, 5]
        job.db.chromatogram.return_value = [1, 2, 3]
        job.db.extractedIonChromatogram.return_value = [1, 2, 3]
        job.db.fragments.return_value = [1, 2, 3]
        return job


class ViewsTestCase(AbstractViewsTestCase):

    def test_home(self):
        request = testing.DummyRequest()
        views = Views(request)
        response = views.home()

        self.assertEqual(response, {})

    def test_allinone_with_ms_data_as_file(self):
        from cgi import FieldStorage
        ms_file = FieldStorage()
        ms_file.filename = r'c:\bla\bla\F1234.mzxml'
        post = {'ms_data_file': ms_file}
        request = testing.DummyRequest(post=post)
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        job = self.fake_job()
        jobquery = Mock(JobQuery)
        job.jobquery.return_value = jobquery
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        views.job_factory.fromScratch = Mock(return_value=job)

        response = views.allinone()

        views.job_factory.fromScratch.assert_called_with('bob')
        jobquery.allinone.assert_called_with(post)
        views.job_factory.submitQuery.assert_called_with(jobquery.allinone(),
                                                         job)
        self.assertEqual(response, {'success': True, 'jobid': 'foo'})
        self.assertEqual(job.ms_filename, r'c:\bla\bla\F1234.mzxml')
        job.jobquery.assert_called_with('http://example.com/status/foo.json',
                                        False,
                                        1)

    def test_allinone_with_ms_data_as_text(self):
        post = {'ms_data': 'somexml', 'ms_data_file': ''}
        request = testing.DummyRequest(post=post)
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        job = self.fake_job()
        jobquery = Mock(JobQuery)
        job.jobquery.return_value = jobquery
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        views.job_factory.fromScratch = Mock(return_value=job)

        response = views.allinone()

        views.job_factory.fromScratch.assert_called_with('bob')
        jobquery.allinone.assert_called_with(post)
        views.job_factory.submitQuery.assert_called_with(jobquery.allinone(),
                                                         job)
        self.assertEqual(response, {'success': True, 'jobid': 'foo'})
        self.assertEqual(job.ms_filename, 'Uploaded as text')
        job.jobquery.assert_called_with('http://example.com/status/foo.json',
                                        False,
                                        1)

    def test_allinone_restricted(self):
        self.config.add_settings(restricted=True)
        post = {'ms_data': 'somexml', 'ms_data_file': ''}
        request = testing.DummyRequest(post=post)
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        job = self.fake_job()
        jobquery = Mock(JobQuery)
        job.jobquery.return_value = jobquery
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        views.job_factory.fromScratch = Mock(return_value=job)

        views.allinone()

        job.jobquery.assert_called_with('http://example.com/status/foo.json',
                                        True,
                                        1)

    def test_allinone_no_jobmanager(self):
        from magmaweb.job import JobSubmissionError
        from pyramid.httpexceptions import HTTPInternalServerError
        import json
        post = {'ms_data': 'somexml', 'ms_data_file': ''}
        request = testing.DummyRequest(post=post)
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        job = self.fake_job()
        jobquery = Mock(JobQuery)
        job.jobquery.return_value = jobquery
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        views.job_factory.fromScratch = Mock(return_value=job)
        q = Mock(side_effect=JobSubmissionError())
        views.job_factory.submitQuery = q

        with self.assertRaises(HTTPInternalServerError) as e:
            views.allinone()

        expected_json = {'success': False, 'msg': 'Unable to submit query'}
        self.assertEqual(json.loads(e.exception.body), expected_json)
        views.job_factory.fromScratch.assert_called_with('bob')
        jobquery.allinone.assert_called_with(post)
        views.job_factory.submitQuery.assert_called_with(jobquery.allinone(),
                                                         job)

    def test_uploaddb_get(self):
        request = testing.DummyRequest()
        views = Views(request)

        response = views.uploaddb()
        self.assertEqual(response, {})

    def test_uploaddb_post(self):
        self.config.add_route('results', '/results/{jobid}')
        from cgi import FieldStorage
        dbfile = FieldStorage()
        dbfile.file = StringIO()
        request = testing.DummyRequest(post={'db_file': dbfile})
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        job = self.fake_job()
        views.job_factory.fromDb.return_value = job

        response = views.uploaddb()

        views.job_factory.fromDb.assert_called_with(dbfile.file, 'bob')
        self.assertEqual(response.location, 'http://example.com/results/foo')

    def test_jobfromscratch(self):
        self.config.add_route('results', '/results/{jobid}')
        request = testing.DummyRequest()
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        job = self.fake_job()
        views = Views(request)
        views.job_factory.fromScratch = Mock(return_value=job)

        response = views.jobfromscratch()

        views.job_factory.fromScratch.assert_called_with('bob')
        self.assertEqual(response.location, 'http://example.com/results/foo')

    def test_failed_validation(self):
        from colander import Invalid
        from pyramid.httpexceptions import HTTPInternalServerError
        import json
        e = Mock(Invalid)
        e.asdict.return_value = {'query': 'Bad query field',
                                 'format': 'Something wrong in form'
                                 }
        request = testing.DummyRequest()
        request.exception = e
        # use alternate view callable argument convention
        # because exception is passed as context
        views = Views(request)

        response = views.failed_validation()

        self.assertIsInstance(response, HTTPInternalServerError)
        expected = {'success': False,
                    'errors': {'query': 'Bad query field',
                               'format': 'Something wrong in form'
                               }
                    }
        self.assertEqual(json.loads(response.body), expected)

    def test_workspace(self):
        import uuid
        self.config.add_route('results', '/results/{jobid}')
        request = testing.DummyRequest()
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        created_at = datetime.datetime(2012, 11, 14, 10, 48, 26, 504478)
        jobs = [JobMeta(uuid.UUID('11111111-1111-1111-1111-111111111111'),
                        'bob', description='My job', created_at=created_at,
                        ms_filename='F1234.mzxml',
                        state='STOPPED',
                        is_public=False,)]
        request.user.jobs = jobs
        views = Views(request)
        views.job_factory.dbSize = Mock(return_value=1234)

        response = views.workspace()

        id1 = '11111111-1111-1111-1111-111111111111'
        url1 = 'http://example.com/results/' + id1
        expected_jobs = [{'id': '11111111-1111-1111-1111-111111111111',
                          'url': url1,
                          'description': 'My job',
                          'is_public': False,
                          'ms_filename': 'F1234.mzxml',
                          'state': 'STOPPED',
                          'created_at': '2012-11-14T10:48:26',
                          'size': 1234,
                          }]
        self.assertEqual(response, {'jobs': expected_jobs})

    @patch('magmaweb.views.JobQuery')
    def test_defaultsjson(self, jq):
        request = testing.DummyRequest()
        views = Views(request)
        jq.defaults.return_value = 'foo'

        response = views.defaults()

        jq.defaults.assert_called_with(None)
        self.assertDictEqual(response, {'success': True, 'data': 'foo'})

    @patch('magmaweb.views.JobQuery')
    def test_examplejson(self, jq):
        request = testing.DummyRequest(params={'selection': 'example'})
        views = Views(request)
        jq.defaults.return_value = 'foo'

        response = views.defaults()

        jq.defaults.assert_called_with('example')
        self.assertDictEqual(response, {'success': True, 'data': 'foo'})

    def test_login_get_from_loginpage(self):
        self.config.add_route('home', '/')
        self.config.add_route('login', '/login')
        request = testing.DummyRequest()
        request.user = None
        request.url = 'http://example.com/login'
        route_mapper = self.config.get_routes_mapper()
        request.matched_route = route_mapper.get_route('login')
        views = Views(request)

        response = views.login()

        expected_response = {'came_from': 'http://example.com/',
                             'userid': '',
                             'password': ''
                             }
        self.assertDictEqual(response, expected_response)
        self.assertEqual(request.response.headers['WWW-Authenticate'], 'MAC')
        self.assertEqual(request.response.status_int, 401)

    def test_login_get(self):
        self.config.add_route('home', '/')
        self.config.add_route('login', '/login')
        self.config.add_route('workspace', '/workspace')
        request = testing.DummyRequest()
        request.user = None
        request.url = 'http://example.com/workspace'
        route_mapper = self.config.get_routes_mapper()
        request.matched_route = route_mapper.get_route('workspace')
        views = Views(request)

        response = views.login()

        expected_response = {'came_from': 'http://example.com/workspace',
                             'userid': '',
                             'password': ''
                             }
        self.assertDictEqual(response, expected_response)
        self.assertEqual(request.response.headers['WWW-Authenticate'], 'MAC')
        self.assertEqual(request.response.status_int, 401)

    @patch('magmaweb.views.remember')
    @patch('magmaweb.views.User')
    def test_login_post_correct_pw(self, user, remember):
        from pyramid.httpexceptions import HTTPFound
        u = Mock(User)
        u.validate_password.return_value = True
        user.by_id.return_value = u
        self.config.add_route('login', '/login')
        self.config.add_route('workspace', '/workspace')
        post = {'userid': 'bob',
                'password': 'mypw'}
        params = {'came_from': 'http://example.com/workspace'}
        request = testing.DummyRequest(post=post, params=params)
        request.user = None
        route_mapper = self.config.get_routes_mapper()
        request.matched_route = route_mapper.get_route('workspace')
        views = Views(request)

        response = views.login()

        self.assertIsInstance(response, HTTPFound)
        self.assertEqual(response.location, 'http://example.com/workspace')
        u.validate_password.assert_called_with('mypw')
        remember.assert_called_with(request, 'bob')

    @patch('magmaweb.views.remember')
    @patch('magmaweb.views.User')
    def test_login_post_incorrect_pw(self, user, remember):
        u = Mock(User)
        u.validate_password.return_value = False
        user.by_id.return_value = u
        self.config.add_route('login', '/login')
        self.config.add_route('workspace', '/workspace')
        post = {'userid': 'bob',
                'password': 'mypw'}
        params = {'came_from': 'http://example.com/workspace'}
        request = testing.DummyRequest(post=post, params=params)
        request.user = None
        route_mapper = self.config.get_routes_mapper()
        request.matched_route = route_mapper.get_route('workspace')

        views = Views(request)

        response = views.login()

        expected_response = {'came_from': 'http://example.com/workspace',
                             'userid': 'bob',
                             'password': 'mypw'
                             }
        self.assertDictEqual(response, expected_response)
        u.validate_password.assert_called_with('mypw')
        self.assertFalse(remember.called)

    def test_login_authenticated(self):
        from pyramid.exceptions import Forbidden
        request = testing.DummyRequest()
        request.user = 'bob'
        request.exception = Forbidden()
        views = Views(request)

        response = views.login()

        self.assertIsInstance(response, Forbidden)

    @patch('magmaweb.views.remember')
    @patch('magmaweb.views.User')
    def test_login_auto_register(self, user, remember):
        user.generate.return_value = User('bob',
                                          'bob dob',
                                          'bob.dob@example.com'
                                          )
        from pyramid.httpexceptions import HTTPFound
        self.config.add_route('home', '/')
        self.config.add_route('login', '/login')
        self.config.add_route('workspace', '/workspace')
        request = testing.DummyRequest()
        request.user = None
        request.url = 'http://example.com/workspace'
        route_mapper = self.config.get_routes_mapper()
        request.matched_route = route_mapper.get_route('workspace')
        request.registry.settings['auto_register'] = True
        views = Views(request)

        response = views.login()

        self.assertIsInstance(response, HTTPFound)
        self.assertEqual(response.location, 'http://example.com/workspace')
        user.generate.assert_called_with()
        remember.assert_called_with(request, 'bob')

    @patch('magmaweb.views.remember')
    @patch('magmaweb.views.User')
    def test_login_auto_register_from_login(self, user, remember):
        user.generate.return_value = User('bob',
                                          'bob dob',
                                          'bob.dob@example.com'
                                          )
        from pyramid.httpexceptions import HTTPFound
        self.config.add_route('home', '/')
        self.config.add_route('login', '/login')
        request = testing.DummyRequest()
        request.user = None
        request.url = 'http://example.com/login'
        route_mapper = self.config.get_routes_mapper()
        request.matched_route = route_mapper.get_route('login')
        request.registry.settings['auto_register'] = True
        views = Views(request)

        response = views.login()

        self.assertIsInstance(response, HTTPFound)
        self.assertEqual(response.location, 'http://example.com/')
        user.generate.assert_called_with()
        remember.assert_called_with(request, 'bob')

    def test_login_auto_register_mac_challenge(self):
        self.config.add_route('status.json', '/status')
        self.config.add_route('home', '/')
        self.config.add_route('login', '/login')
        request = testing.DummyRequest()
        request.user = None
        request.url = 'http://example.com/status'
        request.method = 'PUT'
        route_mapper = self.config.get_routes_mapper()
        request.matched_route = route_mapper.get_route('status.json')
        request.registry.settings['auto_register'] = True
        views = Views(request)

        response = views.login()

        expected_response = {'came_from': 'http://example.com/status',
                             'userid': '',
                             'password': ''
                             }
        self.assertDictEqual(response, expected_response)
        self.assertEqual(request.response.headers['WWW-Authenticate'], 'MAC')
        self.assertEqual(request.response.status_int, 401)

    @patch('magmaweb.views.forget')
    def test_logout(self, forget):
        self.config.add_route('home', '/')
        forget.return_value = {'forget': 'headers'}
        request = testing.DummyRequest()
        views = Views(request)

        response = views.logout()

        forget.assert_called_with(request)
        self.assertEqual(response.location, 'http://example.com/')
        self.assertDictContainsSubset({'forget': 'headers'}, response.headers)

    @patch('magmaweb.views.time.time')
    def test_access_token(self, time):
        time.return_value = 1355000000.000000
        request = testing.DummyRequest()
        request.user = User('bob', 'Bob Example', 'bob@example.com')

        # mock auth
        from pyramid.authorization import ACLAuthorizationPolicy
        from pyramid_macauth import MACAuthenticationPolicy
        from pyramid_multiauth import MultiAuthenticationPolicy
        self.config.set_authorization_policy(ACLAuthorizationPolicy())
        policy = Mock(MACAuthenticationPolicy)
        policy.encode_mac_id.return_value = 'id', 'key'
        mpolicy = MultiAuthenticationPolicy([policy])
        self.config.set_authentication_policy(mpolicy)

        views = Views(request)

        response = views.access_token()

        expected_response = {"acesss_token": 'id',
                             "mac_key": 'key',
                             "expires_in": 360,
                             "token_type": "mac",
                             "mac_algorithm": "hmac-sha-1"}
        self.assertDictEqual(response, expected_response)
        self.assertEqual(request.response.headers['Cache-Control'], "no-store")
        policy.encode_mac_id.assert_called_with(request, 'bob',
                                                expires=1355000360)

    def test_help(self):
        request = testing.DummyRequest()
        views = Views(request)
        response = views.help()

        self.assertEqual(response, {})


class InCompleteJobViewsTestCase(AbstractViewsTestCase):

    """Test case for magmaweb.views.InCompleteJobViews"""

    def test_jobstatus_complete(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'RUNNING'
        job.is_complete.side_effect = JobIncomplete(job)
        views = InCompleteJobViews(job, request)

        response = views.job_status()

        self.assertEqual(response, dict(status='RUNNING',
                                        jobid='bla',
                                        is_complete=False,
                                        ))

    def test_jobstatus_incomplete(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'STOPPED'
        job.is_complete.return_value = True
        views = InCompleteJobViews(job, request)

        response = views.job_status()

        self.assertEqual(response, dict(status='STOPPED',
                                        jobid='bla',
                                        is_complete=True,
                                        ))

    def test_jobstatus_completewitherror(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'ERROR'
        job.is_complete.side_effect = JobError(job)
        views = InCompleteJobViews(job, request)

        with self.assertRaises(JobError) as e:
            views.job_status()

        self.assertEqual(e.exception.job, job)

    def test_jobstatusjson_complete(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'RUNNING'
        job.is_complete.side_effect = JobIncomplete(job)
        views = InCompleteJobViews(job, request)

        response = views.job_status_json()

        self.assertEqual(response, dict(status='RUNNING',
                                        jobid='bla',
                                        is_complete=False,
                                        ))

    def test_jobstatusjson_incomplete(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'STOPPED'
        job.is_complete.return_value = True
        views = InCompleteJobViews(job, request)

        response = views.job_status_json()

        self.assertEqual(response, dict(status='STOPPED',
                                        jobid='bla',
                                        is_complete=True,
                                        ))

    def test_jobstatusjson_completewitherror(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'ERROR'
        job.is_complete.side_effect = JobError(job)
        views = InCompleteJobViews(job, request)

        response = views.job_status_json()

        self.assertEqual(response, dict(status='ERROR',
                                        jobid='bla',
                                        is_complete=True,
                                        ))

    def test_set_jobstatus(self):
        request = testing.DummyRequest(content_type='text/plain')
        request.body = 'STOPPED'
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'RUNNING'
        views = InCompleteJobViews(job, request)

        response = views.set_job_status()

        self.assertEqual(response, dict(status='STOPPED', jobid='bla'))
        self.assertEqual(job.state, 'STOPPED')

    def test_set_jobstatus_withjson_running(self):
        body = '{"state": "EXECUTING", "exitCode": null, "running": true, '
        body += '"done": false, "schedulerSpecficInformation": null, '
        body += '"exception": null}'
        request = testing.DummyRequest(content_type='application/json',
                                       charset='UTF-8',
                                       body=body,
                                       )
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'PENDING'
        views = InCompleteJobViews(job, request)

        response = views.set_job_status()

        self.assertEqual(response, dict(status='EXECUTING', jobid='bla'))
        self.assertEqual(job.state, 'EXECUTING')

    def test_set_jobstatus_withjson_doneok(self):
        body = '{"state": "DONE", "exitCode": 0, "running": false, '
        body += '"done": true, "schedulerSpecficInformation": null, '
        body += '"exception": null}'
        request = testing.DummyRequest(content_type='application/json',
                                       charset='UTF-8',
                                       body=body,
                                       )
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'PENDING'
        views = InCompleteJobViews(job, request)

        response = views.set_job_status()

        self.assertEqual(response, dict(status='STOPPED', jobid='bla'))
        self.assertEqual(job.state, 'STOPPED')

    def test_set_jobstatus_withjson_doneerr_killed(self):
        body = '{"state": "KILLED", "exitCode": 0, "running": false, '
        body += '"done": true, "schedulerSpecficInformation": null, '
        body += '"exception": "local adaptor: Process cancelled by user."}'
        request = testing.DummyRequest(content_type='application/json',
                                       charset='UTF-8',
                                       body=body,
                                       )
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'PENDING'
        views = InCompleteJobViews(job, request)

        response = views.set_job_status()

        self.assertEqual(response, dict(status='ERROR', jobid='bla'))
        self.assertEqual(job.state, 'ERROR')

    def test_set_jobstatus_withjson_doneerr_exitcode(self):
        body = '{"state": "STOPPED", "exitCode": 1, "running": false, '
        body += '"done": true, "schedulerSpecficInformation": null, '
        body += '"exception": null}'
        request = testing.DummyRequest(content_type='application/json',
                                       charset='UTF-8',
                                       body=body,
                                       )
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'PENDING'
        views = InCompleteJobViews(job, request)

        response = views.set_job_status()

        self.assertEqual(response, dict(status='ERROR', jobid='bla'))
        self.assertEqual(job.state, 'ERROR')

    def test_updatejson(self):
        request = testing.DummyRequest()
        request.json_body = {"id": "bar",
                             "description": "New description",
                             "ms_filename": "F12345.mzxml",
                             "created_at": "1999-12-17T13:45:04",
                             "is_public": True,
                             }
        job = self.fake_job()
        expected_id = job.id
        execpted_ca = job.created_at
        views = InCompleteJobViews(job, request)

        response = views.updatejson()

        exp_response = {'success': True, 'message': 'Updated job'}
        self.assertDictEqual(response, exp_response)
        self.assertEqual(job.id, expected_id)
        self.assertEqual(job.description, 'New description')
        self.assertEqual(job.ms_filename, 'F12345.mzxml')
        self.assertEqual(job.created_at, execpted_ca)
        self.assertEqual(job.is_public, True)

    def test_delete_completedjob(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        views = InCompleteJobViews(job, request)

        response = views.delete()

        exp_response = {'success': True, 'message': 'Deleted job'}
        self.assertDictEqual(response, exp_response)
        job.delete.assert_called_with()
        self.assertEqual(request.response.status_int, 204)

    def test_delete_incompletejob(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        # make job incomplete
        job.is_complete.side_effect = JobIncomplete(job)
        views = InCompleteJobViews(job, request)
        # see if cancel is called on job_factory
        views.job_factory = Mock(JobFactory)

        views.delete()

        views.job_factory.cancel.assert_called_with(job)

    def test_delete_erroredjob(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        # make job incomplete
        job.is_complete.side_effect = JobError(job)
        views = InCompleteJobViews(job, request)
        # see if cancel is called on job_factory
        views.job_factory = Mock(JobFactory)

        views.delete()

        job.delete.assert_called_with()

    def test_delete_unable2cancel(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        # make job incomplete
        job.is_complete.side_effect = JobIncomplete(job)
        views = InCompleteJobViews(job, request)
        from requests.exceptions import ConnectionError
        exc = ConnectionError('[Errno 111] Connection refused')
        views.job_factory.cancel = Mock(side_effect=exc)
        from pyramid.httpexceptions import HTTPInternalServerError

        with self.assertRaises(HTTPInternalServerError) as e:
            views.delete()

        expected = json.loads('{"msg": "Failed to cancel job", "success": false}')
        value = json.loads(e.exception.body)
        self.assertEqual(value, expected)

    def test_error(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        exc = JobError(job)
        views = InCompleteJobViews(exc, request)

        response = views.error()

        eresponse = {'exception': exc,
                     'job': job,
                     'run': job.db.runInfo(),
                     }

        self.assertEqual(response, eresponse)

    def test_job_incomplete(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'RUNNING'
        job.is_complete.return_value = False
        exc = JobIncomplete(job)
        views = InCompleteJobViews(exc, request)

        response = views.job_incomplete()

        self.assertEqual(response, dict(status='RUNNING',
                                        jobid='bla',
                                        is_complete=False,
                                        ))

    @patch('magmaweb.views.has_permission')
    def test_resuls_canrun(self, has_permission):
        from pyramid.security import Allowed
        has_permission.return_value = Allowed('Faked allowed')
        request = testing.DummyRequest()
        job = self.fake_job()
        views = InCompleteJobViews(job, request)
        response = views.results()

        self.assertEqual(response, dict(jobid='foo',
                                        run='bla',
                                        maxmslevel=3,
                                        canRun=True,
                                        job=job,
                                        ))
        has_permission.assert_called_with('run', job, request)

    @patch('magmaweb.views.has_permission')
    def test_resuls_cantrun(self, has_permission):
        from pyramid.security import Denied
        has_permission.return_value = Denied('Faked denied')
        request = testing.DummyRequest()
        job = self.fake_job()
        views = InCompleteJobViews(job, request)

        response = views.results()

        exp_response = dict(jobid='foo',
                            run='bla',
                            maxmslevel=3,
                            canRun=False,
                            job=job,
                            )
        self.assertDictEqual(response, exp_response)
        has_permission.assert_called_with('run', job, request)

    def test_stderr(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        log = StringIO()
        log.write('bla')
        log.seek(0)
        job.stderr.return_value = log
        views = InCompleteJobViews(job, request)

        response = views.stderr()

        self.assertEqual(response.content_type, 'text/plain')
        self.assertMultiLineEqual(response.app_iter.read(), 'bla')

    def test_stdout(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        log = StringIO()
        log.write('bla')
        log.seek(0)
        job.stdout.return_value = log
        views = InCompleteJobViews(job, request)

        response = views.stdout()

        self.assertEqual(response.content_type, 'text/plain')
        self.assertMultiLineEqual(response.app_iter.read(), 'bla')


class JobViewsTestCase(AbstractViewsTestCase):

    """ Test case for magmaweb.views.JobViews"""

    def setUp(self):
        AbstractViewsTestCase.setUp(self)
        self.config.add_route('status', '/status/{jobid}')

    def test_incompletejob(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        job = self.fake_job()
        job.is_complete = Mock()

        JobViews(job, request)

        job.is_complete.assert_called_with()

    def test_moleculesjson_return(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        response = views.moleculesjson()

        self.assertEqual(response, {'totalUnfiltered': 3,
                                    'total': 3,
                                    'rows': [1, 2, 3],
                                    'scans': [4, 5]
                                    })

    def test_moleculesjson_minimalparams(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.moleculesjson()

        job.db.molecules.assert_called_with(start=0,
                                            limit=10,
                                            sorts=[],
                                            scanid=None,
                                            filters=[],
                                            molid=None,
                                            mz=None,
                                            )
        job.db.scansWithMolecules.assert_called_with(filters=[], molid=None, mz=None, scanid=None)

    def test_moleculesjson_scanidfilter(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'scanid': 641
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.moleculesjson()

        job.db.molecules.assert_called_with(start=0,
                                            limit=10,
                                            sorts=[],
                                            scanid=641,
                                            filters=[],
                                            molid=None,
                                            mz=None,
                                            )

    def test_moleculesjson_nrscaneqfilter(self):
        filter_in = '[{"type":"numeric","comparison":"eq",'
        filter_in += '"value":1,"field":"nhits"}]'
        filter_expected = [{"type": "numeric", "comparison": "eq",
                            "value": 1, "field": "nhits"}]
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'filter': filter_in
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.moleculesjson()

        job.db.molecules.assert_called_with(start=0,
                                            limit=10,
                                            sorts=[],
                                            scanid=None,
                                            filters=filter_expected,
                                            molid=None,
                                            mz=None,
                                            )
        job.db.scansWithMolecules.assert_called_with(filters=filter_expected,
                                                     molid=None,
                                                     mz=None,
                                                     scanid=None)

    def test_moleculesjson_notfilledreturn(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        job = self.fake_job()
        job.db.moleculesTotalCount.return_value = 0
        views = JobViews(job, request)

        response = views.moleculesjson()

        self.assertEqual(response, {'totalUnfiltered': 0,
                                    'total': 3,
                                    'rows': [1, 2, 3],
                                    'scans': [4, 5]})

    def test_moleculesjson_molpage(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'molid': 1,
                                               })
        job = self.fake_job()
        views = JobViews(job, request)
        job.db.molecules.return_value = {'total': 3, 'rows': [1, 2, 3], 'page': 1}

        response = views.moleculesjson()

        self.assertEqual(response, {'totalUnfiltered': 3,
                                    'total': 3,
                                    'rows': [1, 2, 3],
                                    'molid': 1,
                                    'page': 1,
                                    'scans': [4, 5]})

    def test_moleculesjson_scansfilteredonmol(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'molid': 1,
                                               })
        job = self.fake_job()
        views = JobViews(job, request)
        job.db.molecules.return_value = {'total': 3, 'rows': [1, 2, 3], 'page': 1, 'molid': 1}

        views.moleculesjson()

        job.db.scansWithMolecules.assert_called_once_with(filters=[], molid=1, mz=None, scanid=None)

    def test_moleculescsv(self):
        csv = StringIO()
        csv.write('bla')
        job = self.fake_job()
        job.db.molecules2csv.return_value = csv
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        views = JobViews(job, request)
        views.moleculesjson = Mock(return_value={'rows': []})

        response = views.moleculescsv()

        self.assertEqual(response.content_type, 'text/csv')
        self.assertEqual(response.body, b'bla')

    def test_moleculescsv_somecols(self):
        csv = StringIO()
        csv.write('bla')
        job = self.fake_job()
        job.db.molecules2csv.return_value = csv
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'cols': '["name","score"]'
                                               })

        views = JobViews(job, request)
        rows = [{'name': 'foo', 'score': 'bar', 'id': 123}]
        views.moleculesjson = Mock(return_value={'rows': rows})

        views.moleculescsv()

        job.db.molecules2csv.assert_called_with(rows, cols=['name', 'score'])

    def test_moleculessdf(self):
        job = self.fake_job()
        job.db.molecules2sdf.return_value = 'bla'
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        views = JobViews(job, request)
        views.moleculesjson = Mock(return_value={'rows': []
                                                 })
        response = views.moleculessdf()

        job.db.molecules2sdf.assert_called_with([], cols=[])
        self.assertEqual(response.content_type, 'chemical/x-mdl-sdfile')
        self.assertEqual(response.body, b'bla')

    def test_moleculessdf_somecols(self):
        job = self.fake_job()
        job.db.molecules2sdf.return_value = 'bla'
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'cols': '["name","score"]'
                                               })
        views = JobViews(job, request)
        views.moleculesjson = Mock(return_value={'rows': []
                                                 })
        views.moleculessdf()

        job.db.molecules2sdf.assert_called_with([], cols=['name', 'score'])

    def test_chromatogramjson(self):
        request = testing.DummyRequest()
        views = JobViews(self.fake_job(), request)

        response = views.chromatogramjson()

        self.assertEqual(response, [1, 2, 3])

    def test_mspectrajson(self):
        request = testing.DummyRequest(matchdict={'scanid': 641})
        job = self.fake_job()
        views = JobViews(job, request)

        views.mspectrajson()

        job.db.mspectra.assert_called_with(641, None)

    def test_mspectrajson_withmslevel(self):
        request = testing.DummyRequest(matchdict={'scanid': 641},
                                       params={'mslevel': 3}
                                       )
        job = self.fake_job()
        views = JobViews(job, request)

        views.mspectrajson()

        job.db.mspectra.assert_called_with(641, 3)

    def test_mspectra_withoutscanid_notfound(self):
        from magmaweb.job import ScanNotFound
        request = testing.DummyRequest(matchdict={'scanid': 641})

        job = self.fake_job()
        job.db.mspectra.side_effect = ScanNotFound()
        views = JobViews(job, request)

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            views.mspectrajson()

    def test_extractedionchromatogram(self):
        request = testing.DummyRequest(matchdict={'molid': 72})
        job = self.fake_job()
        views = JobViews(job, request)

        response = views.extractedionchromatogram()

        self.assertEqual(response, {
            'chromatogram': [1, 2, 3], 'scans': [4, 5]
        })
        job.db.extractedIonChromatogram.assert_called_with(72)
        job.db.scansWithMolecules.assert_called_with(molid=72)

    def test_fragments_moleculewithoutfragments(self):
        request = testing.DummyRequest(matchdict={'molid': 72,
                                                  'scanid': 641
                                                  },
                                       params={'node': ''})
        job = self.fake_job()
        views = JobViews(job, request)

        views.fragments()

        job.db.fragments.assert_called_with(molid=72, scanid=641, node='')

    def test_fragments_filteronmslevel(self):
        request = testing.DummyRequest(matchdict={'molid': 72,
                                                  'scanid': 641
                                                  },
                                       params={'node': 1709})
        job = self.fake_job()
        views = JobViews(job, request)

        views.fragments()

        job.db.fragments.assert_called_with(molid=72, scanid=641, node=1709)

    def test_fragments_badmolecule_notfound(self):
        from magmaweb.job import FragmentNotFound
        request = testing.DummyRequest(matchdict={'molid': 72,
                                                  'scanid': 641
                                                  },
                                       params={'node': ''})

        job = self.fake_job()
        job.db.fragments.side_effect = FragmentNotFound
        views = JobViews(job, request)

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            views.fragments()

    def test_runinfojson(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        from magmaweb.models import Run
        job.db.runInfo.return_value = Run(
            ionisation_mode=-1, skip_fragmentation=True,
            ms_intensity_cutoff=200000.0, msms_intensity_cutoff=50,
            mz_precision=4.0, mz_precision_abs=0.002, use_all_peaks=True,
            ms_filename='F123456.mzxml', abs_peak_cutoff=1000,
            max_ms_level=3, precursor_mz_precision=0.01,
            max_broken_bonds=4, description='My first description',
            max_water_losses=1,
        )
        views = JobViews(job, request)

        response = views.runinfojson()

        self.assertEqual(response, {'success': True,
                                    'data': dict(ms_data_area='',
                                                 ms_data_format='mzxml',
                                                 ionisation_mode=-1,
                                                 ms_intensity_cutoff=200000.0,
                                                 msms_intensity_cutoff=50,
                                                 mz_precision=4.0,
                                                 mz_precision_abs=0.002,
                                                 abs_peak_cutoff=1000,
                                                 max_ms_level=3,
                                                 precursor_mz_precision=0.01,
                                                 max_broken_bonds=4,
                                                 max_water_losses=1,
                                                 )
                                    })

    def test_runinfojson_onlymsdatadone(self):
        self.maxDiff = 200000
        request = testing.DummyRequest()
        job = self.fake_job()
        from magmaweb.models import Run
        job.db.runInfo.return_value = Run(abs_peak_cutoff=1100)
        views = JobViews(job, request)

        response = views.runinfojson()

        self.assertEqual(response, {'success': True,
                                    'data': dict(ms_data_area='',
                                                 ms_data_format='mzxml',
                                                 ionisation_mode=1,
                                                 ms_intensity_cutoff=0.0,
                                                 msms_intensity_cutoff=5,
                                                 mz_precision=5.0,
                                                 mz_precision_abs=0.001,
                                                 abs_peak_cutoff=1100,
                                                 max_ms_level=10,
                                                 precursor_mz_precision=0.005,
                                                 max_broken_bonds=3,
                                                 max_water_losses=1,
                                                 )
                                    })

    def test_runinfojson_norundone(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.db.runInfo.return_value = None
        views = JobViews(job, request)

        response = views.runinfojson()

        self.assertEqual(response, {'success': True,
                                    'data': dict(ms_data_area='',
                                                 ms_data_format='mzxml',
                                                 ionisation_mode=1,
                                                 ms_intensity_cutoff=0.0,
                                                 msms_intensity_cutoff=5,
                                                 mz_precision=5.0,
                                                 mz_precision_abs=0.001,
                                                 abs_peak_cutoff=5000,
                                                 max_ms_level=10,
                                                 precursor_mz_precision=0.005,
                                                 max_broken_bonds=3,
                                                 max_water_losses=1,
                                                 )
                                    })
