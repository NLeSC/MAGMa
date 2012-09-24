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
        self.assertEquals(response, { 'success': True, 'jobid': self.jobid})

    def test_addmsdata(self):
        self.jobquery.add_ms_data.return_value = self.jq

        response = self.rpc.add_ms_data()

        self.jobquery.add_ms_data.assert_called_with(self.post, True)
        self.job.metabolitesTotalCount.assert_called_with()
        self.rpc.new_job.assert_called_with()
        self.job.jobquery.assert_called_with()
        self.rpc.job_factory.submitQuery.assert_called_with(self.jq)
        self.assertEquals(response, { 'success': True, 'jobid': self.jobid})

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
        self.assertEquals(response, { 'success': True, 'jobid': self.jobid})

    def test_set_description(self):
        self.rpc.request.POST = { 'description': 'My description'}
        self.rpc.job = Mock(return_value=self.job)

        response = self.rpc.set_description()

        self.job.description.assert_called_with('My description')
        self.assertEquals(response, { 'success': True, 'jobid': self.jobid})

    def test_submit_query_withoutjobmanager(self):
        from urllib2 import URLError
        from pyramid.httpexceptions import HTTPInternalServerError
        import json

        self.rpc.job_factory.submitQuery = Mock(side_effect=URLError('[Errno 111] Connection refused'))
        with self.assertRaises(HTTPInternalServerError) as e:
            self.rpc.submit_query({})

        self.assertEquals(json.loads(e.exception.body), { 'success': False, 'msg': 'Unable to submit query'})

    def test_failed_validation(self):
        from colander import Invalid
        from pyramid.httpexceptions import HTTPInternalServerError
        import json
        e = Mock(Invalid)
        e.asdict.return_value = {
                                 'query': 'Bad query field',
                                 'format': 'Something wrong in form'
                                 }
        request = testing.DummyRequest()
        # use alternate view callable argument convention
        # because exception is passed as context
        rpc = RpcViews(e, request)

        response = rpc.failed_validation()

        self.assertIsInstance(response, HTTPInternalServerError)
        self.assertEqual(json.loads(response.body), {
                                    'success': False,
                                    'errors': {
                                               'query': 'Bad query field',
                                               'format': 'Something wrong in form'
                                               }
                                    })

    def test_assign_metabolite2peak(self):
        self.rpc.request.POST = { 'scanid': 641, 'mz': 109.029563903808, 'metid': 72}
        self.rpc.job = Mock(return_value=self.job)

        response = self.rpc.assign_metabolite2peak()

        self.job.assign_metabolite2peak.assert_called_with(641, 109.029563903808, 72)
        self.assertEquals(response, { 'success': True, 'jobid': self.jobid })

    def test_unassign_metabolite2peak(self):
        self.rpc.request.POST = { 'scanid': 641, 'mz': 109.029563903808}
        self.rpc.job = Mock(return_value=self.job)

        response = self.rpc.unassign_metabolite2peak()

        self.job.unassign_metabolite2peak.assert_called_with(641, 109.029563903808)
        self.assertEquals(response, { 'success': True, 'jobid': self.jobid })
