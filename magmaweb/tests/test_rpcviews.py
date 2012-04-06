import unittest
from pyramid import testing
from mock import Mock
from magmaweb.rpc import RpcViews
from magmaweb.job import JobFactory, Job, JobQuery

class RpcViewsTestCase(unittest.TestCase):
    def setUp(self):
        self.settings = {
                         'jobfactory.root_dir': '/somedir',
                         }
        self.config = testing.setUp(settings=self.settings)

        self.jobid = '3ad25048-26f6-11e1-851e-00012e260790'
        self.post = {'key': 'value'}
        self.rpc = RpcViews(testing.DummyRequest(post=self.post))
        self.jq = Mock(JobQuery)
        self.jobquery = Mock(JobQuery)
        self.job = Mock(Job)
        self.job.id = self.jobid
        self.job.jobquery.return_value = self.jobquery
        self.job.maxMSLevel.return_value = 1
        self.job.metabolitesTotalCount.return_value = 1
        self.rpc.job = self.job
        self.rpc.new_job = Mock(return_value=self.job)
        self.rpc.job_factory.submitQuery = Mock()

    def tearDown(self):
        testing.tearDown()

    def test_construct(self):
        request = testing.DummyRequest()

        rpc = RpcViews(request)

        self.assertIsInstance(rpc.job_factory, JobFactory)
        self.assertEqual(request, rpc.request)

    def test_job(self):
        jobid = '3ad25048-26f6-11e1-851e-00012e260790'
        request = testing.DummyRequest()
        request.matchdict['jobid'] = jobid
        rpc = RpcViews(request)
        rpc.job_factory.fromId = Mock(return_value=12345)

        job = rpc.job()

        self.assertEqual(job, 12345)
        rpc.job_factory.fromId.assert_called_with(jobid)

    def test_new_job_with_jobid(self):
        jobid = '3ad25048-26f6-11e1-851e-00012e260790'
        request = testing.DummyRequest()
        request.matchdict['jobid'] = jobid
        rpc = RpcViews(request)
        rpc.job_factory.fromId = Mock(return_value=12345)
        rpc.job_factory.cloneJob = Mock(return_value=67890)

        job = rpc.new_job()

        self.assertEqual(job, 67890)
        rpc.job_factory.fromId.assert_called_with(jobid)
        rpc.job_factory.cloneJob.assert_called_with(12345)

    def test_new_job_without_jobid(self):
        request = testing.DummyRequest()
        rpc = RpcViews(request)
        rpc.job_factory.fromScratch = Mock(return_value=12345)

        job = rpc.new_job()

        self.assertEqual(job, 12345)
        rpc.job_factory.fromScratch.assert_called_with()

    def test_addstructure(self):
        self.jobquery.add_structures.return_value = self.jq

        response = self.rpc.add_structures()

        self.jobquery.add_structures.assert_called_with(self.post, True)
        self.job.maxMSLevel.assert_called_with()
        self.rpc.new_job.assert_called_with()
        self.job.jobquery.assert_called_with()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEqual(response.content_type, 'text/html') # required for extjs iframe file upload
        self.assertEqual(response.body, '{"success": true, "jobid": "'+str(self.jobid)+'"}')

    def test_addmsdata(self):
        self.jobquery.add_ms_data.return_value = self.jq

        response = self.rpc.add_ms_data()

        self.jobquery.add_ms_data.assert_called_with(self.post, True)
        self.job.metabolitesTotalCount.assert_called_with()
        self.rpc.new_job.assert_called_with()
        self.job.jobquery.assert_called_with()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEqual(response.content_type, 'text/html') # required for extjs iframe file upload
        self.assertEqual(response.body, '{"success": true, "jobid": "'+str(self.jobid)+'"}')

    def test_metabolize(self):
        self.jobquery.metabolize.return_value = self.jq

        response = self.rpc.metabolize()

        self.jobquery.metabolize.assert_called_with(self.post, True)
        self.job.maxMSLevel.assert_called_with()
        self.rpc.new_job.assert_called_with()
        self.job.jobquery.assert_called_with()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEquals(response, { 'success': True, 'jobid': self.jobid})

    def test_metabolize_one(self):
        self.jobquery.metabolize_one.return_value = self.jq

        response = self.rpc.metabolize_one()

        self.jobquery.metabolize_one.assert_called_with(self.post, True)
        self.job.maxMSLevel.assert_called_with()
        self.rpc.new_job.assert_called_with()
        self.job.jobquery.assert_called_with()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEquals(response, { 'success': True, 'jobid': self.jobid})

    def test_annotate(self):
        self.jobquery.annotate.return_value = self.jq

        response = self.rpc.annotate()

        self.jobquery.annotate.assert_called_with(self.post)
        self.rpc.new_job.assert_called_with()
        self.job.jobquery.assert_called_with()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEquals(response, { 'success': True, 'jobid': self.jobid})

    def test_allinone(self):
        self.jobquery.allinone.return_value = self.jq

        response = self.rpc.allinone()

        self.jobquery.allinone.assert_called_with(self.post)
        self.rpc.new_job.assert_called_with()
        self.job.jobquery.assert_called_with()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEqual(response.content_type, 'text/html') # required for extjs iframe file upload
        self.assertEqual(response.body, '{"success": true, "jobid": "'+str(self.jobid)+'"}')
