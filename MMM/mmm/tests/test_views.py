import unittest
from pyramid import testing
from mock import patch

class HelperTestCase(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.settings = { 'jobrootdir': tempfile.mkdtemp() }
        self.config = testing.setUp(settings=self.settings)
        self.config.add_route('results', '/results')
        self.config.add_route('home', '/homepage')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.config.registry.settings['jobrootdir'])
        testing.tearDown()

    def test_job_factory(self):
        from mmm.views import job_factory
        request = testing.DummyRequest()
        jf = job_factory(request)
        self.assertEqual(jf.jobrootdir, self.settings['jobrootdir'])
        self.assertEqual(jf.dbname, 'results.db')

    @patch('mmm.views.job_factory')
    def test_fetch_job_session(self, mocked_jobfactory):
        from mmm.job import JobFactory
        from mock import Mock
        jobf = Mock(JobFactory)
        job = mock_job()
        job.id = 'foobar'
        jobf.fromId.return_value = job
        mocked_jobfactory.return_value = jobf
        request = testing.DummyRequest()
        request.session['id'] = job.id

        from mmm.views import fetch_job
        out = fetch_job(request)

        self.assertEqual(out, job)
        jobf.fromId.assert_called_with(request.session['id'])

    @patch('mmm.views.job_factory')
    def test_fetch_job_with_jobid_in_get(self, mocked_jobfactory):
        from mmm.job import JobFactory
        from mock import Mock
        jobf = Mock(JobFactory)
        job = mock_job()
        job.id = 'foobar'
        jobf.fromId.return_value = job
        mocked_jobfactory.return_value = jobf
        request = testing.DummyRequest()
        request.params['jobid'] = job.id

        from mmm.views import fetch_job
        out = fetch_job(request)

        self.assertEqual(out, job)
        jobf.fromId.assert_called_with(request.session['id'])

    def test_fetch_job_without_id(self):
        request = testing.DummyRequest()

        from mmm.views import fetch_job
        from pyramid.httpexceptions import HTTPFound
        with self.assertRaises(HTTPFound):
            fetch_job(request)

class HomeView(unittest.TestCase):
    def setUp(self):
        import tempfile
        settings = { 'jobrootdir': tempfile.mkdtemp() }
        self.config = testing.setUp(settings=settings)
        self.config.add_route('results', '/results')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.config.registry.settings['jobrootdir'])
        testing.tearDown()

    def test_get(self):
        request = testing.DummyRequest()
        from mmm.views import home
        response = home(request)
        self.assertEqual(response, {})

    @patch('mmm.views.job_factory')
    def test_post(self, mocked_jobfactory):
        job = mock_job()
        job.id = 'x'
        from mmm.job import JobFactory
        from mock import Mock
        jobf = Mock(JobFactory)
        jobf.fromQuery.return_value = job
        mocked_jobfactory.return_value = jobf

        import os
        dbfile = os.tmpfile()

        class FileUpload:
            pass
        post = { 'db': FileUpload() }
        post['db'].file = dbfile
        request = testing.DummyRequest(post=post)

        from mmm.views import home, job_factory
        response = home(request)

        self.assertEqual(response.body, '{"success": true}')
        self.assertEqual(response.content_type, 'text/html')
        self.assertEqual(request.session['id'], job.id)
        mocked_jobfactory.assert_called_with(request)
        jobf.fromQuery.assert_called_with(dbfile)

def mock_job():
    """ Returns a mocked mmm.job.Job which can be used as return value in a patched fetch_job()

    Example:
    @patch('mmm.views.fetch_job')
    def test_it(self, mocked_fetch_job):
        job = mock_job()
        job.maxMSLevel.return_value = 3
        mocked_fetch_job.return_value = job

        ... do action

        assert job.maxMSLevel.called

    """
    from mmm.job import Job
    from mock import Mock
    job = Mock(Job)
    return job

class ResultsView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()

    @patch('mmm.views.fetch_job')
    def test_it(self, mocked_fetch_job):
        job = mock_job()
        job.runInfo.return_value = 'bla'
        job.maxMSLevel.return_value = 3
        mocked_fetch_job.return_value = job

        request = testing.DummyRequest()
        from mmm.views import results, fetch_job
        response = results(request)

        self.assertEqual(response, dict(run='bla', maxmslevel=3))
        job.runInfo.assert_called_with()
        job.maxMSLevel.assert_called_with()

class MetabolitesView(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()
        self.job = mock_job()
        self.job.metabolites.return_value = {'total':3, 'rows':[1,2,3]}
        self.job.scansWithMetabolites.return_value = [4,5]

    def tearDown(self):
        testing.tearDown()

    def _callFUT(self, params):
        from mmm.views import metabolitesjson, fetch_job
        request = testing.DummyRequest(params=params)
        return metabolitesjson(request)

    @patch('mmm.views.fetch_job')
    def test_default(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        response = self._callFUT(params=dict(start=0, limit=10))

        self.assertEqual(response, { 'total':3, 'rows':[1,2,3], 'scans':[4,5]})
        self.job.metabolites.assert_called_with(start=0, limit=10, sorts=[],scanid=None, filters=[])
        self.job.scansWithMetabolites.assert_called_with(filters=[])

    @patch('mmm.views.fetch_job')
    def test_filteredon_scanid(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        params = dict(start=0, limit=10, scanid=641)
        response = self._callFUT(params)

        self.assertEqual(response, { 'total':3, 'rows':[1,2,3], 'scans':[4,5]})
        self.job.metabolites.assert_called_with(start=0, limit=10, scanid=641, sorts=[], filters=[])
        self.job.scansWithMetabolites.assert_called_with(filters=[])

    @patch('mmm.views.fetch_job')
    def test_filteredon_nrscanseq(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        params = dict(start=0, limit=10, filter='[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}]')
        response = self._callFUT(params)

        self.assertEqual(response, { 'total':3, 'rows':[1,2,3], 'scans':[4,5]})
        self.job.metabolites.assert_called_with(start=0, limit=10, scanid=None, sorts=[], filters=[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}])
        self.job.scansWithMetabolites.assert_called_with(filters=[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}])

    @patch('mmm.views.fetch_job')
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

    @patch('mmm.views.fetch_job')
    def test_it(self,  mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        from mmm.views import metabolitescsv, fetch_job
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
        from mmm.views import chromatogramjson, fetch_job
        request = testing.DummyRequest(params=params)
        return chromatogramjson(request)

    @patch('mmm.views.fetch_job')
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
        from mmm.views import mspectrajson
        request = testing.DummyRequest(params=params)
        request.matchdict = {'scanid': id }
        return mspectrajson(request)

    @patch('mmm.views.fetch_job')
    def test_scanid(self, mocked_fetch_job):
        self.job.mspectra.return_value = 'foobar'
        mocked_fetch_job.return_value = self.job

        self.assertEqual(self._callFUT(641,dict()), 'foobar')
        self.job.mspectra.assert_called_with(641, None)

    @patch('mmm.views.fetch_job')
    def test_scanid_and_mslevel(self, mocked_fetch_job):
        self.job.mspectra.return_value = 'foobar'
        mocked_fetch_job.return_value = self.job

        self.assertEqual(self._callFUT(641,dict(mslevel=1)), 'foobar')
        self.job.mspectra.assert_called_with(641, 1)

    @patch('mmm.views.fetch_job')
    def test_notfound(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job
        from mmm.job import ScanNotFound
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
        from mmm.views import extractedionchromatogram
        request = testing.DummyRequest()
        request.matchdict = {'metid': id }
        return extractedionchromatogram(request)

    @patch('mmm.views.fetch_job')
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
        from mmm.views import fragments
        request = testing.DummyRequest(params=params)
        request.matchdict = {'metid': metid, 'scanid': scanid }
        return fragments(request)

    @patch('mmm.views.fetch_job')
    def test_metabolitewithoutfragments(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        response = self._callFUT(72, 641, dict(node=''))

        self.assertEqual(response, [1, 2, 3])
        self.job.fragments.assert_called_with(metid=72, scanid=641, node='')

    @patch('mmm.views.fetch_job')
    def test_lvl3fragments(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job

        response = self._callFUT(352, 870, dict(node=1709))

        self.assertEqual(response, [1, 2, 3])
        self.job.fragments.assert_called_with(metid=352, scanid=870, node=1709)

    @patch('mmm.views.fetch_job')
    def test_badmetabolite(self, mocked_fetch_job):
        mocked_fetch_job.return_value = self.job
        from mmm.job import FragmentNotFound
        self.job.fragments.side_effect = FragmentNotFound()

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            self._callFUT(70002, 641, dict(node=''))
