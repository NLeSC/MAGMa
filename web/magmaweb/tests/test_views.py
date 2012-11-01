import unittest
from pyramid import testing
from mock import Mock, patch
from magmaweb.views import Views, JobViews
from magmaweb.job import JobFactory, Job, JobDb

class AbstractViewsTestCase(unittest.TestCase):
    def setUp(self):
        self.settings = {
                         'jobfactory.root_dir': '/somedir',
                         }
        self.config = testing.setUp(settings=self.settings)

    def tearDown(self):
        testing.tearDown()

    def fake_job(self):
        job = Mock(Job)
        job.id = 'foo'
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

        self.assertEqual(response, {'userid':None})

    def test_uploaddb_get(self):
        request = testing.DummyRequest()
        views = Views(request)

        response = views.uploaddb()
        self.assertEqual(response, {})

    @patch('magmaweb.views.unauthenticated_userid')
    def test_uploaddb_post(self, unauthenticated_userid):
        unauthenticated_userid.return_value = 'bob'
        self.config.add_route('results', '/results/{jobid}')
        from cgi import FieldStorage
        import StringIO
        dbfile = FieldStorage()
        dbfile.file = StringIO.StringIO()
        request = testing.DummyRequest(post={'db_file': dbfile})
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        job = self.fake_job()
        views.job_factory.fromDb.return_value = job

        response = views.uploaddb()

        views.job_factory.fromDb.assert_called_with(dbfile.file, 'bob')
        self.assertEqual(response.location, 'http://example.com/results/foo')

    @patch('magmaweb.views.unauthenticated_userid')
    def test_jobfromscratch(self, unauthenticated_userid):
        unauthenticated_userid.return_value = 'bob'
        self.config.add_route('results', '/results/{jobid}')
        request = testing.DummyRequest()
        job = self.fake_job()
        views = Views(request)
        views.job_factory.fromScratch = Mock(return_value=job)

        response = views.jobfromscratch()

        views.job_factory.fromScratch.assert_called_with('bob')
        self.assertEqual(response.location, 'http://example.com/results/foo')

    def test_results(self):
        request = testing.DummyRequest()
        views = JobViews(self.fake_job(),request)

        response = views.results()

        self.assertEqual(response, dict(
                                        jobid='foo',
                                        run='bla',
                                        maxmslevel=3
                                        ))

    @patch('magmaweb.views.get_jobs')
    def test_jobs(self, get_jobs):
        jobs = [{'id':12345, 'description': 'My job'}]
        get_jobs.return_value = jobs
        request = testing.DummyRequest()
        views = Views(request)

        response = views.jobs()

        get_jobs.assert_called_with(None)
        self.assertEqual(response, {'jobs':jobs})

    def test_defaultsjson(self):
        request = testing.DummyRequest()
        views = Views(request)
        response = views.defaults()

        self.assertEqual(response, {
                                    'success': True,
                                    'data': dict(
                                                 n_reaction_steps=2,
                                                 metabolism_types=['phase1',
                                                                   'phase2'],
                                                 ionisation_mode=1,
                                                 skip_fragmentation=False,
                                                 ms_intensity_cutoff=1000000.0,
                                                 msms_intensity_cutoff=0.1,
                                                 mz_precision=0.001,
                                                 use_all_peaks=False,
                                                 abs_peak_cutoff=1000,
                                                 rel_peak_cutoff=0.01,
                                                 max_ms_level=10,
                                                 precursor_mz_precision=0.005,
                                                 max_broken_bonds=4
                                                 )
                                    })


class JobViewsTestCase(AbstractViewsTestCase):
    """ Test case for magmaweb.views.JobViews"""
    def test_jobstatus(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.id = 'bla'
        job.state = 'RUNNING'
        views = JobViews(job, request)

        response = views.job_status()

        self.assertEqual(response, dict(status='RUNNING', jobid='bla'))

    def test_metabolitesjson_return(self):
        request = testing.DummyRequest(params={
                                               'start': 0,
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
        request = testing.DummyRequest(params={
                                               'start': 0,
                                               'limit': 10
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.metabolitesjson()

        job.db.metabolites.assert_called_with(start=0, limit=10, sorts=[],
                                           scanid=None, filters=[])
        job.db.scansWithMetabolites.assert_called_with(filters=[])

    def test_metabolitesjson_scanidfilter(self):
        request = testing.DummyRequest(params={
                                               'start': 0,
                                               'limit': 10,
                                               'scanid': 641
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.metabolitesjson()

        job.db.metabolites.assert_called_with(start=0, limit=10,
                                           sorts=[], scanid=641,
                                           filters=[])

    def test_metabolitesjson_nrscaneqfilter(self):
        filter_in = '[{"type":"numeric","comparison":"eq",'
        filter_in += '"value":1,"field":"nhits"}]'
        filter_expected = [{"type": "numeric", "comparison": "eq",
                            "value": 1, "field": "nhits"}]
        request = testing.DummyRequest(params={
                                               'start': 0,
                                               'limit': 10,
                                               'filter': filter_in
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.metabolitesjson()

        job.db.metabolites.assert_called_with(start=0, limit=10,
                                           sorts=[], scanid=None,
                                           filters=filter_expected)
        job.db.scansWithMetabolites.assert_called_with(filters=filter_expected)

    def test_metabolitesjson_sortonlevel(self):
        sort_in = '[{"property":"level","direction":"DESC"}]'
        sort_expected = [{"property":"level", "direction": "DESC"}]
        request = testing.DummyRequest(params={
                                               'start': 0,
                                               'limit': 10,
                                               'sort': sort_in
                                               })
        job = self.fake_job()
        views = JobViews(job, request)

        views.metabolitesjson()

        job.db.metabolites.assert_called_with(start=0, limit=10,
                                           sorts=sort_expected,
                                           scanid=None, filters=[])

    def test_metabolitesjson_notfilledreturn(self):
        request = testing.DummyRequest(params={
                                               'start': 0,
                                               'limit': 10
                                               })
        job = self.fake_job()
        job.db.metabolitesTotalCount.return_value = 0
        views = JobViews(job, request)

        response = views.metabolitesjson()

        self.assertEqual(response, {'totalUnfiltered': 0, 'total': 3,
                                    'rows': [1, 2, 3], 'scans': [4, 5]})

    def test_metabolitescsv(self):
        import StringIO
        csv = StringIO.StringIO()
        csv.write('bla')
        job = self.fake_job()
        job.db.metabolites2csv.return_value = csv
        request = testing.DummyRequest(params={
                                               'start': 0,
                                               'limit': 10
                                               })
        views = JobViews(job, request)
        views.metabolitesjson = Mock(return_value={
                                                   'rows': []
                                                   })

        response = views.metabolitescsv()

        self.assertEqual(response.content_type, 'text/csv')
        self.assertEqual(response.body, 'bla')

    def test_metabolitescsv_somecols(self):
        import StringIO
        csv = StringIO.StringIO()
        csv.write('bla')
        job = self.fake_job()
        job.db.metabolites2csv.return_value = csv
        request = testing.DummyRequest(params={
                                               'start': 0,
                                               'limit': 10,
                                               'cols': '["name","score"]'
                                               })

        views = JobViews(job, request)
        rows = [{'name':'foo', 'score': 'bar', 'id': 123}]
        views.metabolitesjson = Mock(return_value={'rows': rows})

        views.metabolitescsv()

        job.db.metabolites2csv.assert_called_with(rows, cols=['name', 'score'])

    def test_metabolitessdf(self):
        job = self.fake_job()
        job.db.metabolites2sdf.return_value = 'bla'
        request = testing.DummyRequest(params={
                                               'start': 0,
                                               'limit': 10
                                               })
        views = JobViews(job, request)
        views.metabolitesjson = Mock(return_value={
                                                   'rows': []
                                                   })
        response = views.metabolitessdf()

        job.db.metabolites2sdf.assert_called_with([], cols=[])
        self.assertEqual(response.content_type, 'chemical/x-mdl-sdfile')
        self.assertEqual(response.body, 'bla')

    def test_metabolitessdf_somecols(self):
        job = self.fake_job()
        job.db.metabolites2sdf.return_value = 'bla'
        request = testing.DummyRequest(params={
                                               'start': 0,
                                               'limit': 10,
                                               'cols': '["name","score"]'
                                               })
        views = JobViews(job, request)
        views.metabolitesjson = Mock(return_value={
                                                   'rows': []
                                                   })
        views.metabolitessdf()

        job.db.metabolites2sdf.assert_called_with([], cols=['name', 'score'])

    def test_chromatogramjson(self):
        request = testing.DummyRequest()
        views = JobViews(self.fake_job(), request)

        response = views.chromatogramjson()

        self.assertEqual(response, [1, 2, 3])

    def test_mspectrajson(self):
        request = testing.DummyRequest(matchdict={'scanid':641})
        job = self.fake_job()
        views = JobViews(job, request)

        views.mspectrajson()

        job.db.mspectra.assert_called_with(641, None)

    def test_mspectrajson_withmslevel(self):
        request = testing.DummyRequest(
                                       matchdict={'scanid': 641},
                                       params={'mslevel': 3}
                                       )
        job = self.fake_job()
        views = JobViews(job, request)

        views.mspectrajson()

        job.db.mspectra.assert_called_with(641, 3)

    def test_mspectra_withoutscanid_notfound(self):
        from magmaweb.job import ScanNotFound
        request = testing.DummyRequest(matchdict={'scanid':641})

        job = self.fake_job()
        job.db.mspectra.side_effect = ScanNotFound()
        views = JobViews(job, request)

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            views.mspectrajson()

    def test_extractedionchromatogram(self):
        request = testing.DummyRequest(matchdict={'metid':72})
        job = self.fake_job()
        views = JobViews(job, request)

        response = views.extractedionchromatogram()

        self.assertEqual(response, {
            'chromatogram': [1, 2, 3], 'scans': [4, 5]
        })
        job.db.extractedIonChromatogram.assert_called_with(72)
        job.db.scansWithMetabolites.assert_called_with(metid=72)

    def test_fragments_metabolitewithoutfragments(self):
        request = testing.DummyRequest(matchdict={
                                                  'metid': 72,
                                                  'scanid': 641
                                                  },
                                       params={'node':''})
        job = self.fake_job()
        views = JobViews(job, request)

        views.fragments()

        job.db.fragments.assert_called_with(metid=72, scanid=641, node='')

    def test_fragments_filteronmslevel(self):
        request = testing.DummyRequest(matchdict={
                                                  'metid': 72,
                                                  'scanid': 641
                                                  },
                                       params={'node':1709})
        job = self.fake_job()
        views = JobViews(job, request)

        views.fragments()

        job.db.fragments.assert_called_with(metid=72, scanid=641, node=1709)

    def test_fragments_badmetabolite_notfound(self):
        from magmaweb.job import FragmentNotFound
        request = testing.DummyRequest(matchdict={
                                                  'metid': 72,
                                                  'scanid': 641
                                                  },
                                       params={'node':''})

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
            n_reaction_steps=2, metabolism_types='phase1,phase2',
            ionisation_mode=-1, skip_fragmentation=True,
            ms_intensity_cutoff=200000.0, msms_intensity_cutoff=0.5,
            mz_precision=0.01, use_all_peaks=True,
            ms_filename='F123456.mzxml', abs_peak_cutoff=1000,
            rel_peak_cutoff=0.001, max_ms_level=3, precursor_mz_precision=0.01,
            max_broken_bonds=4, description='My first description'
        )
        views = JobViews(job, request)

        response = views.runinfojson()

        self.assertEqual(response, {
                                    'success': True,
                                    'data': dict(
                                                 n_reaction_steps=2,
                                                 metabolism_types=['phase1',
                                                                   'phase2'],
                                                 ionisation_mode=-1,
                                                 skip_fragmentation=True,
                                                 ms_intensity_cutoff=200000.0,
                                                 msms_intensity_cutoff=0.5,
                                                 mz_precision=0.01,
                                                 use_all_peaks=True,
                                                 abs_peak_cutoff=1000,
                                                 rel_peak_cutoff=0.001,
                                                 max_ms_level=3,
                                                 precursor_mz_precision=0.01,
                                                 max_broken_bonds=4
                                                 )
                                    })

    def test_runinfojson_onlymsdatadone(self):
        self.maxDiff = 200000
        request = testing.DummyRequest()
        job = self.fake_job()
        from magmaweb.models import Run
        job.db.runInfo.return_value = Run(
            abs_peak_cutoff=1100,
            rel_peak_cutoff=0.012
        )
        views = JobViews(job, request)

        response = views.runinfojson()

        self.assertEqual(response, {
                                    'success': True,
                                    'data': dict(
                                                 n_reaction_steps=2,
                                                 metabolism_types=['phase1',
                                                                   'phase2'],
                                                 ionisation_mode=1,
                                                 skip_fragmentation=False,
                                                 ms_intensity_cutoff=1000000.0,
                                                 msms_intensity_cutoff=0.1,
                                                 mz_precision=0.001,
                                                 use_all_peaks=False,
                                                 abs_peak_cutoff=1100,
                                                 rel_peak_cutoff=0.012,
                                                 max_ms_level=10,
                                                 precursor_mz_precision=0.005,
                                                 max_broken_bonds=4
                                                 )
                                    })

    def test_runinfojson_norundone(self):
        request = testing.DummyRequest()
        job = self.fake_job()
        job.db.runInfo.return_value = None
        views = JobViews(job, request)

        response = views.runinfojson()

        self.assertEqual(response, {
                                    'success': True,
                                    'data': dict(
                                                 n_reaction_steps=2,
                                                 metabolism_types=['phase1',
                                                                   'phase2'],
                                                 ionisation_mode=1,
                                                 skip_fragmentation=False,
                                                 ms_intensity_cutoff=1000000.0,
                                                 msms_intensity_cutoff=0.1,
                                                 mz_precision=0.001,
                                                 use_all_peaks=False,
                                                 abs_peak_cutoff=1000,
                                                 rel_peak_cutoff=0.01,
                                                 max_ms_level=10,
                                                 precursor_mz_precision=0.005,
                                                 max_broken_bonds=4
                                                 )
                                    })
