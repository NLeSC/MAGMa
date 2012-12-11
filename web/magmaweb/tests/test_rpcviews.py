import unittest
import uuid
from pyramid import testing
from mock import Mock
from magmaweb.rpc import RpcViews
from magmaweb.job import JobFactory, Job, JobDb, JobQuery
from magmaweb.user import JobMeta, User


class RpcViewsTestCase(unittest.TestCase):
    def setUp(self):
        self.settings = {'jobfactory.root_dir': '/somedir',
                         }
        self.config = testing.setUp(settings=self.settings)
        self.config.add_route('status.json', '/status/{jobid}.json')

        self.jobid = '3ad25048-26f6-11e1-851e-00012e260790'
        self.post = {'key': 'value'}
        self.jq = Mock(JobQuery)
        self.jobquery = Mock(JobQuery)
        jobmeta = JobMeta(uuid.UUID(self.jobid), owner='bob')
        self.job = Job(jobmeta, '/mydir', Mock(JobDb))
        self.job.jobquery = Mock(return_value=self.jobquery)
        self.job.db.maxMSLevel.return_value = 1
        self.job.db.metabolitesTotalCount.return_value = 1

        self.rpc = RpcViews(self.job, testing.DummyRequest(post=self.post))
        self.rpc.new_job = Mock(return_value=self.job)
        self.rpc.job_factory.submitQuery = Mock()

    def tearDown(self):
        testing.tearDown()

    def _assert_status_callback_url(self):
        expected_cb_url = 'http://example.com/status/'
        expected_cb_url += '3ad25048-26f6-11e1-851e-00012e260790.json'
        self.job.jobquery.assert_called_with(expected_cb_url)

    def test_construct(self):
        request = testing.DummyRequest()

        rpc = RpcViews(self.job, request)

        self.assertIsInstance(rpc.job_factory, JobFactory)
        self.assertEqual(request, rpc.request)
        self.assertEqual(self.job, rpc.job)

    def test_new_job(self):
        request = testing.DummyRequest()
        request.user = User('bob', 'Bob Example', 'bob@example.com')
        parent_job = Mock(Job)
        rpc = RpcViews(parent_job, request)
        rpc.job_factory.cloneJob = Mock(return_value=self.job)

        job = rpc.new_job()

        self.assertEqual(job, self.job)
        rpc.job_factory.cloneJob.assert_called_with(parent_job, 'bob')

    def test_addstructure(self):
        self.jobquery.add_structures.return_value = self.jq

        response = self.rpc.add_structures()

        self.jobquery.add_structures.assert_called_with(self.post, True)
        self.job.db.maxMSLevel.assert_called_with()
        self.rpc.new_job.assert_called_with()
        self._assert_status_callback_url()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEquals(response, {'success': True, 'jobid': self.jobid})

    def test_addmsdata(self):
        from cgi import FieldStorage
        ms_file = FieldStorage()
        ms_file.filename = 'c:\bla\bla\F1234.mzxml'
        post = {'ms_data_file': ms_file}
        self.rpc.request.POST = post
        self.jobquery.add_ms_data.return_value = self.jq

        response = self.rpc.add_ms_data()

        self.jobquery.add_ms_data.assert_called_with(post, True)
        self.job.db.metabolitesTotalCount.assert_called_with()
        self.assertEqual(self.job.ms_filename, 'c:\bla\bla\F1234.mzxml')
        self.rpc.new_job.assert_called_with()
        self._assert_status_callback_url()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEquals(response, {'success': True, 'jobid': self.jobid})

    def test_metabolize(self):
        self.jobquery.metabolize.return_value = self.jq

        response = self.rpc.metabolize()

        self.jobquery.metabolize.assert_called_with(self.post, True)
        self.job.db.maxMSLevel.assert_called_with()
        self.rpc.new_job.assert_called_with()
        self._assert_status_callback_url()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEquals(response, {'success': True, 'jobid': self.jobid})

    def test_metabolize_one(self):
        self.jobquery.metabolize_one.return_value = self.jq

        response = self.rpc.metabolize_one()

        self.jobquery.metabolize_one.assert_called_with(self.post, True)
        self.job.db.maxMSLevel.assert_called_with()
        self.rpc.new_job.assert_called_with()
        self._assert_status_callback_url()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEquals(response, {'success': True, 'jobid': self.jobid})

    def test_annotate(self):
        self.jobquery.annotate.return_value = self.jq

        response = self.rpc.annotate()

        self.jobquery.annotate.assert_called_with(self.post)
        self.rpc.new_job.assert_called_with()
        self._assert_status_callback_url()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEquals(response, {'success': True, 'jobid': self.jobid})

    def test_set_description(self):
        self.rpc.request.POST = {'description': 'My description'}

        response = self.rpc.set_description()

        self.assertEquals(self.job.description, 'My description')
        self.assertEquals(response, {'success': True, 'jobid': self.jobid})

    def test_submit_query_withoutjobmanager(self):
        from magmaweb.job import JobSubmissionError
        from pyramid.httpexceptions import HTTPInternalServerError
        import json

        q = Mock(side_effect=JobSubmissionError())
        self.rpc.job_factory.submitQuery = q
        with self.assertRaises(HTTPInternalServerError) as e:
            self.rpc.submit_query({})

        expected_json = {'success': False, 'msg': 'Unable to submit query'}
        self.assertEquals(json.loads(e.exception.body), expected_json)

    def test_assign_metabolite2peak(self):
        self.rpc.request.POST = {'scanid': 641,
                                 'mz': 109.029563903808,
                                 'metid': 72}

        response = self.rpc.assign_metabolite2peak()

        assigned_func = self.job.db.assign_metabolite2peak
        assigned_func.assert_called_with(641, 109.029563903808, 72)
        self.assertEquals(response, {'success': True, 'jobid': self.jobid})

    def test_unassign_metabolite2peak(self):
        self.rpc.request.POST = {'scanid': 641,
                                 'mz': 109.029563903808}

        response = self.rpc.unassign_metabolite2peak()

        unassign_func = self.job.db.unassign_metabolite2peak
        unassign_func.assert_called_with(641, 109.029563903808)
        self.assertEquals(response, {'success': True, 'jobid': self.jobid})
