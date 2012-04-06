import unittest
from pyramid import testing
from mock import patch, Mock

class FileUpload:
    file = None
    filename = None
    pass

class HelperTestCase(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.settings = {
                         'jobfactory.root_dir': tempfile.mkdtemp(),
                         }
        self.config = testing.setUp(settings=self.settings)
        self.config.add_route('results', '/results/{jobid}')
        self.config.add_route('home', '/homepage')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.config.registry.settings['jobfactory.root_dir'])
        testing.tearDown()

    def test_job_factory(self):
        from magmaweb.views import job_factory
        request = testing.DummyRequest()
        jf = job_factory(request)
        self.assertEqual(jf.root_dir, self.settings['jobfactory.root_dir'])

    @patch('magmaweb.views.job_factory')
    def test_fetch_job(self, mocked_jobfactory):
        from magmaweb.job import JobFactory
        jobf = Mock(JobFactory)
        job = mock_job()
        job.id = 'foobar'
        jobf.fromId.return_value = job
        mocked_jobfactory.return_value = jobf
        request = testing.DummyRequest()
        request.matchdict['jobid'] = job.id

        from magmaweb.views import fetch_job
        out = fetch_job(request)

        self.assertEqual(out, job)
        jobf.fromId.assert_called_with(request.matchdict['jobid'])

    def test_fetch_job_without_id(self):
        request = testing.DummyRequest()

        from magmaweb.views import fetch_job
        from pyramid.httpexceptions import HTTPFound
        with self.assertRaises(HTTPFound):
            fetch_job(request)

    def test_fetch_job_withwrong_id(self):
        request = testing.DummyRequest()
        request.matchdict['jobid'] = '11111111-1111-1111-1111-111111111111'

        from magmaweb.views import fetch_job
        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            fetch_job(request)

class HomeView(unittest.TestCase):
    def setUp(self):
        import tempfile
        settings = { 'root_dir': tempfile.mkdtemp() }
        self.config = testing.setUp(settings=settings)
        self.config.add_route('results', '/results/{jobid}')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.config.registry.settings['root_dir'])
        testing.tearDown()

    def test_get(self):
        request = testing.DummyRequest()
        from magmaweb.views import home
        response = home(request)
        self.assertEqual(response, {})

    @patch('magmaweb.views.job_factory')
    def test_post(self, mocked_jobfactory):
        import uuid
        jobid = uuid.UUID('3ad25048-26f6-11e1-851e-00012e260790')
        from magmaweb.job import JobFactory, Job, JobQuery
        jobquery = Mock(JobQuery)
        job = Mock(Job)
        job.id = jobid
        job.jobquery.return_value = jobquery
        jobf = Mock(JobFactory)
        jobf.fromScratch.return_value = job
        mocked_jobfactory.return_value = jobf

        post={ 'key': 'value' }
        request = testing.DummyRequest(post=post)

        from magmaweb.views import home
        response = home(request)

        self.assertEqual(response.body, '{"success": true, "jobid": "'+str(jobid)+'"}')
        self.assertEqual(response.content_type, 'text/html') # required for extjs iframe file upload
        mocked_jobfactory.assert_called_with(request)
        jobf.fromScratch.assert_called_with()
        job.jobquery.assert_called_with()
        jobquery.allinone.assert_called_with(post)

def mock_job():
    """ Returns a mocked magmaweb.job.Job which can be used as return value in a patched fetch_job()

    Example:
    @patch('magmaweb.views.fetch_job')
    def test_it(self, mocked_fetch_job):
        job = mock_job()
        job.maxMSLevel.return_value = 3
        mocked_fetch_job.return_value = job

        ... do action

        assert job.maxMSLevel.called

    """
    from magmaweb.job import Job
    job = Mock(Job)
    return job

class ResultsView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()

    @patch('magmaweb.views.fetch_job')
    def test_it(self, mocked_fetch_job):
        job = mock_job()
        job.id = 'foo'
        job.runInfo.return_value = 'bla'
        job.maxMSLevel.return_value = 3
        mocked_fetch_job.return_value = job

        request = testing.DummyRequest()
        request.matchdict['jobid'] = job.id
        from magmaweb.views import results, fetch_job
        response = results(request)

        self.assertEqual(response, dict(jobid=job.id, run='bla', maxmslevel=3))
        job.runInfo.assert_called_with()
        job.maxMSLevel.assert_called_with()

class StatusView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()

    @patch('magmaweb.views.job_factory')
    def test_page(self, jf):
        from magmaweb.job import JobFactory
        mjf = Mock(JobFactory)
        mjf.state = Mock(return_value='STOPPED')
        jf.return_value = mjf
        request = testing.DummyRequest()
        jobid = '3ad25048-26f6-11e1-851e-00012e260790'
        request.matchdict['jobid'] = jobid

        from magmaweb.views import job_status
        response = job_status(request)

        self.assertEqual(response, dict(status='STOPPED', jobid=jobid))
        mjf.state.assert_called_once_with(jobid)

    @patch('magmaweb.views.job_factory')
    def test_json(self, jf):
        from magmaweb.job import JobFactory
        mjf = Mock(JobFactory)
        mjf.state = Mock(return_value='STOPPED')
        jf.return_value = mjf
        request = testing.DummyRequest()
        jobid = '3ad25048-26f6-11e1-851e-00012e260790'
        request.matchdict['jobid'] = jobid

        from magmaweb.views import job_statusjson
        response = job_statusjson(request)

        self.assertEqual(response, dict(status='STOPPED', jobid=jobid))
        mjf.state.assert_called_once_with(jobid)


class MetabolitesView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.job = mock_job()
        self.job.metabolites.return_value = {'total':3, 'rows':[1,2,3]}
        self.job.scansWithMetabolites.return_value = [4,5]

    def tearDown(self):
        testing.tearDown()

    def _callFUT(self, params):
        from magmaweb.views import metabolitesjson, fetch_job
        request = testing.DummyRequest(params=params)
        return metabolitesjson(request)

    @patch('magmaweb.views.fetch_job')
    def test_default(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        response = self._callFUT(params=dict(start=0, limit=10))

        self.assertEqual(response, { 'total':3, 'rows':[1,2,3], 'scans':[4,5]})
        self.job.metabolites.assert_called_with(start=0, limit=10, sorts=[],scanid=None, filters=[])
        self.job.scansWithMetabolites.assert_called_with(filters=[])

    @patch('magmaweb.views.fetch_job')
    def test_filteredon_scanid(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        params = dict(start=0, limit=10, scanid=641)
        response = self._callFUT(params)

        self.assertEqual(response, { 'total':3, 'rows':[1,2,3], 'scans':[4,5]})
        self.job.metabolites.assert_called_with(start=0, limit=10, scanid=641, sorts=[], filters=[])
        self.job.scansWithMetabolites.assert_called_with(filters=[])

    @patch('magmaweb.views.fetch_job')
    def test_filteredon_nrscanseq(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        params = dict(start=0, limit=10, filter='[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}]')
        response = self._callFUT(params)

        self.assertEqual(response, { 'total':3, 'rows':[1,2,3], 'scans':[4,5]})
        self.job.metabolites.assert_called_with(start=0, limit=10, scanid=None, sorts=[], filters=[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}])
        self.job.scansWithMetabolites.assert_called_with(filters=[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}])

    @patch('magmaweb.views.fetch_job')
    def test_sort_score(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        params = dict(start=0, limit=10, scanid=641, sort='[{"property":"score","direction":"DESC"}]')
        response = self._callFUT(params)

        self.assertEqual(response, { 'total':3, 'rows':[1,2,3], 'scans':[4,5]})
        self.job.metabolites.assert_called_with(start=0, limit=10, scanid=641, sorts=[{"property":"score","direction":"DESC"}], filters=[])
        self.job.scansWithMetabolites.assert_called_with(filters=[])

class Metabolites2CsvView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.job = mock_job()
        self.job.metabolites.return_value = {'total':3, 'rows':[{'x':'y'}]}
        import StringIO
        csv = StringIO.StringIO()
        csv.write('bla')
        self.job.metabolites2csv.return_value = csv

    def tearDown(self):
        testing.tearDown()

    @patch('magmaweb.views.fetch_job')
    def test_it(self,  mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        from magmaweb.views import metabolitescsv, fetch_job
        params = dict(start=0, limit=10)
        request = testing.DummyRequest(params=params)
        response = metabolitescsv(request)

        self.assertMultiLineEqual(response.body, 'bla')
        self.assertEqual(response.content_type, 'text/csv')
        self.job.metabolites.assert_called_with(start=0, sorts=[], limit=10, scanid=None, filters=[])
        self.job.metabolites2csv.assert_called_with(self.job.metabolites.return_value['rows'])

class ChromatogramView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.job = mock_job()

    def tearDown(self):
        testing.tearDown()

    def _callFUT(self, params):
        from magmaweb.views import chromatogramjson, fetch_job
        request = testing.DummyRequest(params=params)
        return chromatogramjson(request)

    @patch('magmaweb.views.fetch_job')
    def test_it(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job
        self.job.chromatogram.return_value = [1, 2, 3]

        self.assertEqual(self._callFUT(dict()), [1, 2, 3])
        self.job.chromatogram.assert_called_with()

class MSpectraView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_route('mspectra.json', '/mspectra/{scanid}.json')
        self.job = mock_job()

    def tearDown(self):
        testing.tearDown()

    def _callFUT(self, id, params):
        from magmaweb.views import mspectrajson
        request = testing.DummyRequest(params=params)
        request.matchdict = {'scanid': id }
        return mspectrajson(request)

    @patch('magmaweb.views.fetch_job')
    def test_scanid(self, mocked_fetch_job):
        self.job.mspectra.return_value = 'foobar'
        mocked_fetch_job.return_value = self.job

        self.assertEqual(self._callFUT(641,dict()), 'foobar')
        self.job.mspectra.assert_called_with(641, None)

    @patch('magmaweb.views.fetch_job')
    def test_scanid_and_mslevel(self, mocked_fetch_job):
        self.job.mspectra.return_value = 'foobar'
        mocked_fetch_job.return_value = self.job

        self.assertEqual(self._callFUT(641,dict(mslevel=1)), 'foobar')
        self.job.mspectra.assert_called_with(641, 1)

    @patch('magmaweb.views.fetch_job')
    def test_notfound(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job
        from magmaweb.job import ScanNotFound
        self.job.mspectra.side_effect = ScanNotFound()

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            self._callFUT(123, dict())

class MetaboliteScansView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_route('extractedionchromatogram.json','/extractedionchromatogram/{metid}.json')
        self.job = mock_job()
        self.job.extractedIonChromatogram.return_value = [1,2,3]
        self.job.scansWithMetabolites.return_value = [4,5]

    def tearDown(self):
        testing.tearDown()

    def _callFUT(self, id):
        from magmaweb.views import extractedionchromatogram
        request = testing.DummyRequest()
        request.matchdict = {'metid': id }
        return extractedionchromatogram(request)

    @patch('magmaweb.views.fetch_job')
    def test_it(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        response = self._callFUT(72)

        self.assertEqual(response, {
            'chromatogram': [1,2,3], 'scans': [4,5]
        })
        self.job.extractedIonChromatogram.assert_called_with(72)
        self.job.scansWithMetabolites.assert_called_with(metid=72)

class FragmentsView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_route('fragments.json', '/fragments/{scanid}/{metid}.json')
        self.job = mock_job()
        self.job.fragments.return_value = [1, 2, 3]

    def tearDown(self):
        testing.tearDown()

    def _callFUT(self, metid, scanid, params):
        from magmaweb.views import fragments
        request = testing.DummyRequest(params=params)
        request.matchdict = {'metid': metid, 'scanid': scanid }
        return fragments(request)

    @patch('magmaweb.views.fetch_job')
    def test_metabolitewithoutfragments(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        response = self._callFUT(72, 641, dict(node=''))

        self.assertEqual(response, [1, 2, 3])
        self.job.fragments.assert_called_with(metid=72, scanid=641, node='')

    @patch('magmaweb.views.fetch_job')
    def test_lvl3fragments(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        response = self._callFUT(352, 870, dict(node=1709))

        self.assertEqual(response, [1, 2, 3])
        self.job.fragments.assert_called_with(metid=352, scanid=870, node=1709)

    @patch('magmaweb.views.fetch_job')
    def test_badmetabolite(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job
        from magmaweb.job import FragmentNotFound
        self.job.fragments.side_effect = FragmentNotFound()

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            self._callFUT(70002, 641, dict(node=''))


class StderrView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.job = mock_job()
        import StringIO
        log = StringIO.StringIO()
        log.write('bla')
        log.seek(0)
        self.job.stderr.return_value = log

    def tearDown(self):
        testing.tearDown()

    @patch('magmaweb.views.fetch_job')
    def test_it(self,  mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        from magmaweb.views import stderr, fetch_job
        request = testing.DummyRequest()
        response = stderr(request)

        self.job.stderr.assert_called_with()
        self.assertEqual(response.content_type, 'text/plain')
        self.assertMultiLineEqual(response.app_iter.read(), 'bla')

class UploadDbView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.config.add_route('results', '/results/{jobid}')

    def tearDown(self):
        testing.tearDown()

    def test_get(self):
        request = testing.DummyRequest()
        from magmaweb.views import uploaddb
        response = uploaddb(request)
        self.assertEqual(response, {})

    @patch('magmaweb.views.job_factory')
    def test_post(self, mocked_jobfactory):
        import tempfile
        db_file = tempfile.NamedTemporaryFile()
        post = { 'db_file': FileUpload() }
        post['db_file'].file = db_file
        request = testing.DummyRequest(post=post)
        job = FileUpload()
        job.id = '123456'
        from magmaweb.job import JobFactory
        jobf = Mock(JobFactory)
        jobf.fromDb.return_value = job
        mocked_jobfactory.return_value = jobf

        from magmaweb.views import uploaddb
        response = uploaddb(request)

        from pyramid.httpexceptions import HTTPFound
        self.assertIsInstance(response, HTTPFound)
        self.assertEquals(response.location, 'http://example.com/results/123456')
        mocked_jobfactory.assert_called_with(request)
        jobf.fromDb.assert_called_with(db_file)

        db_file.close()
