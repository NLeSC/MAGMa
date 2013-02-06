import unittest
import datetime
from pyramid import testing
from mock import Mock, patch
from magmaweb.views import Views, JobViews
from magmaweb.job import JobFactory, Job, JobDb, JobQuery
from magmaweb.user import User, JobMeta


class AbstractViewsTestCase(unittest.TestCase):
    def setUp(self):
        self.settings = {'jobfactory.root_dir': '/somedir',
                         'access_token.expires_in': 360,
                         }
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
        job.db = Mock(JobDb)
        job.db.runInfo.return_value = 'bla'
        job.db.maxMSLevel.return_value = 3
        job.db.metabolites.return_value = {'total': 3, 'rows': [1, 2, 3]}
        job.db.metabolitesTotalCount.return_value = 3
        job.db.scansWithMetabolites.return_value = [4, 5]
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

    def test_startjob(self):
        request = testing.DummyRequest()
        views = Views(request)
        response = views.startjob()

        self.assertEqual(response, {})

    def test_allinone_with_ms_data_as_file(self):
        from cgi import FieldStorage
        ms_file = FieldStorage()
        ms_file.filename = 'c:\bla\bla\F1234.mzxml'
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
        self.assertEqual(job.ms_filename, 'c:\bla\bla\F1234.mzxml')
        job.jobquery.assert_called_with('http://example.com/status/foo.json')

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
        job.jobquery.assert_called_with('http://example.com/status/foo.json')

    def test_allinone_with_structure_database(self):
        post = {'ms_data': 'somexml',
                'ms_data_file': '',
                'structure_database': 'pubchem'
                }
        request = testing.DummyRequest(post=post)
        s = request.registry.settings
        s['structure_database.pubchem'] = 'data/pubchem.db'
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        job = self.fake_job()
        jobquery = Mock(JobQuery)
        job.jobquery.return_value = jobquery
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        views.job_factory.fromScratch = Mock(return_value=job)

        response = views.allinone()

        views.job_factory.fromScratch.assert_called_with('bob')
        jobquery.allinone.assert_called_with(post, 'data/pubchem.db')
        views.job_factory.submitQuery.assert_called_with(jobquery.allinone(),
                                                         job)
        self.assertEqual(response, {'success': True, 'jobid': 'foo'})
        self.assertEqual(job.ms_filename, 'Uploaded as text')
        job.jobquery.assert_called_with('http://example.com/status/foo.json')

    def test_allinone_with_empty_structure_database(self):
        post = {'ms_data': 'somexml',
                'ms_data_file': '',
                'structure_database': ''
                }
        request = testing.DummyRequest(post=post)
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        job = self.fake_job()
        jobquery = Mock(JobQuery)
        job.jobquery.return_value = jobquery
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        views.job_factory.fromScratch = Mock(return_value=job)

        views.allinone()

        jobquery.allinone.assert_called_with(post)

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
        self.assertEquals(json.loads(e.exception.body), expected_json)
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
        import StringIO
        dbfile = FieldStorage()
        dbfile.file = StringIO.StringIO()
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
                        ms_filename='F1234.mzxml')]
        request.user.jobs = jobs
        views = Views(request)

        response = views.workspace()

        id1 = '11111111-1111-1111-1111-111111111111'
        url1 = 'http://example.com/results/' + id1
        expected_jobs = [{'id': '11111111-1111-1111-1111-111111111111',
                          'url': url1,
                          'description': 'My job',
                          'ms_filename': 'F1234.mzxml',
                          'created_at': '2012-11-14T10:48:26'}]
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
        request = testing.DummyRequest()
        request.user = None
        request.url = 'http://example.com/startjob'
        views = Views(request)

        response = views.login()

        expected_response = {'came_from': 'http://example.com/startjob',
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
        post = {'userid': 'bob',
                'password': 'mypw'}
        params = {'came_from': 'http://example.com/startjob'}
        request = testing.DummyRequest(post=post, params=params)
        request.user = None
        views = Views(request)

        response = views.login()

        self.assertIsInstance(response, HTTPFound)
        self.assertEqual(response.location, 'http://example.com/startjob')
        u.validate_password.assert_called_with('mypw')
        remember.assert_called_with(request, 'bob')

    @patch('magmaweb.views.remember')
    @patch('magmaweb.views.User')
    def test_login_post_incorrect_pw(self, user, remember):
        u = Mock(User)
        u.validate_password.return_value = False
        user.by_id.return_value = u
        self.config.add_route('login', '/login')
        post = {'userid': 'bob',
                'password': 'mypw'}
        params = {'came_from': 'http://example.com/startjob'}
        request = testing.DummyRequest(post=post, params=params)
        request.user = None
        views = Views(request)

        response = views.login()

        expected_response = {'came_from': 'http://example.com/startjob',
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


class JobViewsTestCase(AbstractViewsTestCase):
    """ Test case for magmaweb.views.JobViews"""

    @patch('magmaweb.views.has_permission')
    def test_resuls_canrun(self, has_permission):
        from pyramid.security import Allowed
        has_permission.return_value = Allowed('Faked allowed')
        request = testing.DummyRequest()
        job = self.fake_job()
        views = JobViews(job, request)

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
        views = JobViews(job, request)

        response = views.results()

        exp_response = dict(jobid='foo',
                            run='bla',
                            maxmslevel=3,
                            canRun=False,
                            job=job,
                            )
        self.assertDictEqual(response, exp_response)
        has_permission.assert_called_with('run', job, request)

    def test_jobstatus(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'RUNNING'
        views = JobViews(job, request)

        response = views.job_status()

        self.assertEqual(response, dict(status='RUNNING', jobid='bla'))

    def test_set_jobstatus(self):
        request = testing.DummyRequest()
        request.body = 'STOPPED'
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'RUNNING'
        views = JobViews(job, request)

        response = views.set_job_status()

        self.assertEqual(response, dict(status='STOPPED', jobid='bla'))
        self.assertEqual(job.state, 'STOPPED')

    def test_metabolitesjson_return(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        response = views.metabolitesjson()

        self.assertEqual(response, {'totalUnfiltered': 3,
                                    'total': 3,
                                    'rows': [1, 2, 3],
                                    'scans': [4, 5]
                                    })

    def test_metabolitesjson_minimalparams(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.metabolitesjson()

        job.db.metabolites.assert_called_with(start=0,
                                              limit=10,
                                              sorts=[],
                                              scanid=None,
                                              filters=[])
        job.db.scansWithMetabolites.assert_called_with(filters=[])

    def test_metabolitesjson_scanidfilter(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'scanid': 641
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.metabolitesjson()

        job.db.metabolites.assert_called_with(start=0,
                                              limit=10,
                                              sorts=[],
                                              scanid=641,
                                              filters=[])

    def test_metabolitesjson_nrscaneqfilter(self):
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

        views.metabolitesjson()

        job.db.metabolites.assert_called_with(start=0,
                                              limit=10,
                                              sorts=[],
                                              scanid=None,
                                              filters=filter_expected)
        job.db.scansWithMetabolites.assert_called_with(filters=filter_expected)

    def test_metabolitesjson_sortonlevel(self):
        sort_in = '[{"property":"level","direction":"DESC"}]'
        sort_expected = [{"property": "level", "direction": "DESC"}]
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'sort': sort_in
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.metabolitesjson()

        job.db.metabolites.assert_called_with(start=0,
                                              limit=10,
                                              sorts=sort_expected,
                                              scanid=None,
                                              filters=[])

    def test_metabolitesjson_notfilledreturn(self):
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        job = self.fake_job()
        job.db.metabolitesTotalCount.return_value = 0
        views = JobViews(job, request)

        response = views.metabolitesjson()

        self.assertEqual(response, {'totalUnfiltered': 0,
                                    'total': 3,
                                    'rows': [1, 2, 3],
                                    'scans': [4, 5]})

    def test_metabolitescsv(self):
        import StringIO
        csv = StringIO.StringIO()
        csv.write('bla')
        job = self.fake_job()
        job.db.metabolites2csv.return_value = csv
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        views = JobViews(job, request)
        views.metabolitesjson = Mock(return_value={'rows': []})

        response = views.metabolitescsv()

        self.assertEqual(response.content_type, 'text/csv')
        self.assertEqual(response.body, 'bla')

    def test_metabolitescsv_somecols(self):
        import StringIO
        csv = StringIO.StringIO()
        csv.write('bla')
        job = self.fake_job()
        job.db.metabolites2csv.return_value = csv
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'cols': '["name","score"]'
                                               })

        views = JobViews(job, request)
        rows = [{'name': 'foo', 'score': 'bar', 'id': 123}]
        views.metabolitesjson = Mock(return_value={'rows': rows})

        views.metabolitescsv()

        job.db.metabolites2csv.assert_called_with(rows, cols=['name', 'score'])

    def test_metabolitessdf(self):
        job = self.fake_job()
        job.db.metabolites2sdf.return_value = 'bla'
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10
                                               })
        views = JobViews(job, request)
        views.metabolitesjson = Mock(return_value={'rows': []
                                                   })
        response = views.metabolitessdf()

        job.db.metabolites2sdf.assert_called_with([], cols=[])
        self.assertEqual(response.content_type, 'chemical/x-mdl-sdfile')
        self.assertEqual(response.body, 'bla')

    def test_metabolitessdf_somecols(self):
        job = self.fake_job()
        job.db.metabolites2sdf.return_value = 'bla'
        request = testing.DummyRequest(params={'start': 0,
                                               'limit': 10,
                                               'cols': '["name","score"]'
                                               })
        views = JobViews(job, request)
        views.metabolitesjson = Mock(return_value={'rows': []
                                                   })
        views.metabolitessdf()

        job.db.metabolites2sdf.assert_called_with([], cols=['name', 'score'])

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
        request = testing.DummyRequest(matchdict={'metid': 72})
        job = self.fake_job()
        views = JobViews(job, request)

        response = views.extractedionchromatogram()

        self.assertEqual(response, {
            'chromatogram': [1, 2, 3], 'scans': [4, 5]
        })
        job.db.extractedIonChromatogram.assert_called_with(72)
        job.db.scansWithMetabolites.assert_called_with(metid=72)

    def test_fragments_metabolitewithoutfragments(self):
        request = testing.DummyRequest(matchdict={'metid': 72,
                                                  'scanid': 641
                                                  },
                                       params={'node': ''})
        job = self.fake_job()
        views = JobViews(job, request)

        views.fragments()

        job.db.fragments.assert_called_with(metid=72, scanid=641, node='')

    def test_fragments_filteronmslevel(self):
        request = testing.DummyRequest(matchdict={'metid': 72,
                                                  'scanid': 641
                                                  },
                                       params={'node': 1709})
        job = self.fake_job()
        views = JobViews(job, request)

        views.fragments()

        job.db.fragments.assert_called_with(metid=72, scanid=641, node=1709)

    def test_fragments_badmetabolite_notfound(self):
        from magmaweb.job import FragmentNotFound
        request = testing.DummyRequest(matchdict={'metid': 72,
                                                  'scanid': 641
                                                  },
                                       params={'node': ''})

        job = self.fake_job()
        job.db.fragments.side_effect = FragmentNotFound
        views = JobViews(job, request)

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            views.fragments()

    def test_stderr(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        import StringIO
        log = StringIO.StringIO()
        log.write('bla')
        log.seek(0)
        job.stderr.return_value = log
        views = JobViews(job, request)

        response = views.stderr()

        self.assertEqual(response.content_type, 'text/plain')
        self.assertMultiLineEqual(response.app_iter.read(), 'bla')

    def test_runinfojson(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        from magmaweb.models import Run
        job.db.runInfo.return_value = Run(
            n_reaction_steps=2, metabolism_types=['phase1', 'phase2'],
            ionisation_mode=-1, skip_fragmentation=True,
            ms_intensity_cutoff=200000.0, msms_intensity_cutoff=0.5,
            mz_precision=4.0, mz_precision_abs=0.002, use_all_peaks=True,
            ms_filename='F123456.mzxml', abs_peak_cutoff=1000,
            max_ms_level=3, precursor_mz_precision=0.01,
            max_broken_bonds=4, description='My first description',
            max_water_losses=1,
        )
        views = JobViews(job, request)

        response = views.runinfojson()

        self.assertEqual(response, {'success': True,
                                    'data': dict(n_reaction_steps=2,
                                                 metabolism_types=['phase1',
                                                                   'phase2'],
                                                 ionisation_mode=-1,
                                                 ms_intensity_cutoff=200000.0,
                                                 msms_intensity_cutoff=0.5,
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
                                    'data': dict(n_reaction_steps=2,
                                                 metabolism_types=['phase1',
                                                                   'phase2'],
                                                 ionisation_mode=1,
                                                 ms_intensity_cutoff=1000000.0,
                                                 msms_intensity_cutoff=0.1,
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
                                    'data': dict(n_reaction_steps=2,
                                                 metabolism_types=['phase1',
                                                                   'phase2'],
                                                 ionisation_mode=1,
                                                 ms_intensity_cutoff=1000000.0,
                                                 msms_intensity_cutoff=0.1,
                                                 mz_precision=5.0,
                                                 mz_precision_abs=0.001,
                                                 abs_peak_cutoff=1000,
                                                 max_ms_level=10,
                                                 precursor_mz_precision=0.005,
                                                 max_broken_bonds=3,
                                                 max_water_losses=1,
                                                 )
                                    })

    def test_updatejson(self):
        request = testing.DummyRequest()
        request.json_body = {"id": "bar",
                             "description": "New description",
                             "ms_filename": "F12345.mzxml",
                             "created_at": "1999-12-17T13:45:04",
                             }
        job = self.fake_job()
        expected_id = job.id
        execpted_ca = job.created_at
        views = JobViews(job, request)

        response = views.updatejson()

        exp_response = {'success': True, 'message': 'Updated job'}
        self.assertDictEqual(response, exp_response)
        self.assertEqual(job.id, expected_id)
        self.assertEqual(job.description, 'New description')
        self.assertEqual(job.ms_filename, 'F12345.mzxml')
        self.assertEqual(job.created_at, execpted_ca)

    def test_deletejson(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        views = JobViews(job, request)

        response = views.deletejson()

        exp_response = {'success': True, 'message': 'Deleted job'}
        self.assertDictEqual(response, exp_response)
        job.delete.assert_called_with()
        self.assertEquals(request.response.status_int, 204)
