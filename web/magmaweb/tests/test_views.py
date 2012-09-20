import unittest
from pyramid import testing
from mock import Mock
from magmaweb.views import Views
from magmaweb.job import JobFactory, JobNotFound, Job

class ViewsTestCase(unittest.TestCase):
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
        job.runInfo.return_value = 'bla'
        job.maxMSLevel.return_value = 3
        job.metabolites.return_value = {'total':3, 'rows':[1,2,3]}
        job.metabolitesTotalCount.return_value = 3
        job.scansWithMetabolites.return_value = [4,5]
        job.chromatogram.return_value = [1,2,3]
        job.extractedIonChromatogram.return_value = [1,2,3]
        job.fragments.return_value = [1, 2, 3]
        return job

    def test_jobid(self):
        jobid = '11111111-1111-1111-1111-111111111111'
        request = testing.DummyRequest()
        request.matchdict['jobid'] = jobid
        views = Views(request)

        actual_jobid = views.jobid()
        self.assertEqual(actual_jobid, jobid)

    def test_jobid_withOutJobidInRequest_redirects2Homepage(self):
        self.config.add_route('home', '/homepage')
        request = testing.DummyRequest()
        views = Views(request)

        from pyramid.httpexceptions import HTTPFound
        with self.assertRaises(HTTPFound) as e:
            views.jobid()
            self.assertEqual(e.location, '/homepage')

    def test_job(self):
        jobid = '11111111-1111-1111-1111-111111111111'
        request = testing.DummyRequest()
        views = Views(request)
        views.jobid = Mock(return_value=jobid)
        views.job_factory.fromId = Mock(return_value=Mock(Job))

        job = views.job()

        self.assertIsInstance(job, Job)

    def test_job_withWrongJobid_notFound(self):
        jobid = '11111111-1111-1111-1111-111111111111'
        request = testing.DummyRequest()
        request.matchdict['jobid'] = jobid
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        views.job_factory.fromId.side_effect = JobNotFound(jobid)

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            views.job()

    def test_home(self):
        request = testing.DummyRequest()
        views = Views(request)

        response = views.home()

        self.assertEqual(response, {})

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
        views = Views(request)
        views.job_factory = Mock(JobFactory)
        views.job_factory.fromDb.return_value = self.fake_job()

        response = views.uploaddb()

        views.job_factory.fromDb.assert_called_with(dbfile.file)
        self.assertEqual(response.location, 'http://example.com/results/foo')

    def test_jobfromscratch(self):
        self.config.add_route('results', '/results/{jobid}')
        request = testing.DummyRequest()
        views = Views(request)
        views.job_factory.fromScratch = Mock(return_value=self.fake_job())

        response = views.jobfromscratch()

        self.assertEqual(response.location, 'http://example.com/results/foo')

    def test_results(self):
        request = testing.DummyRequest()
        views = Views(request)
        views.job = Mock(return_value=self.fake_job())

        response = views.results()

        self.assertEqual(response, dict(
                                        jobid='foo',
                                        run='bla',
                                        maxmslevel=3
                                        ))

    def test_jobstatus(self):
        request = testing.DummyRequest()
        views = Views(request)
        views.jobid = Mock(return_value='bla')
        views.job_factory.state = Mock(return_value='RUNNING')

        response = views.job_status()

        self.assertEqual(response, dict(status='RUNNING', jobid='bla'))

    def test_metabolitesjson_return(self):
        request = testing.DummyRequest(params={
                                               'start':0,
                                               'limit':10
                                               })
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        response = views.metabolitesjson()

        self.assertEqual(response, { 'totalUnfiltered': 3, 'total':3, 'rows':[1,2,3], 'scans':[4,5]})

    def test_metabolitesjson_minimalparams(self):
        request = testing.DummyRequest(params={
                                               'start':0,
                                               'limit':10
                                               })
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        views.metabolitesjson()

        job.metabolites.assert_called_with(start=0, limit=10, sorts=[],scanid=None, filters=[])
        job.scansWithMetabolites.assert_called_with(filters=[])

    def test_metabolitesjson_scanidfilter(self):
        request = testing.DummyRequest(params={
                                               'start':0,
                                               'limit':10,
                                               'scanid': 641
                                               })
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        views.metabolitesjson()

        job.metabolites.assert_called_with(start=0, limit=10, sorts=[],scanid=641, filters=[])

    def test_metabolitesjson_nrscaneqfilter(self):
        filter_in = '[{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}]'
        filter_expected = [{"type":"numeric","comparison":"eq","value":1,"field":"nr_scans"}]
        request = testing.DummyRequest(params={
                                               'start':0,
                                               'limit':10,
                                               'filter': filter_in
                                               })
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        views.metabolitesjson()

        job.metabolites.assert_called_with(start=0, limit=10, sorts=[],scanid=None, filters=filter_expected)
        job.scansWithMetabolites.assert_called_with(filters=filter_expected)

    def test_metabolitesjson_sortonlevel(self):
        sort_in = '[{"property":"level","direction":"DESC"}]'
        sort_expected = [{"property":"level","direction":"DESC"}]
        request = testing.DummyRequest(params={
                                               'start':0,
                                               'limit':10,
                                               'sort': sort_in
                                               })
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        views.metabolitesjson()

        job.metabolites.assert_called_with(start=0, limit=10, sorts=sort_expected,scanid=None, filters=[])

    def test_metabolitesjson_notfilledreturn(self):
        request = testing.DummyRequest(params={
                                               'start':0,
                                               'limit':10
                                               })
        views = Views(request)
        job = self.fake_job()
        job.metabolitesTotalCount.return_value = 0
        views.job = Mock(return_value=job)

        response = views.metabolitesjson()

        self.assertEqual(response, { 'totalUnfiltered': 0, 'total':3, 'rows':[1,2,3], 'scans':[4,5]})

    def test_metabolitescsv(self):
        import StringIO
        csv = StringIO.StringIO()
        csv.write('bla')
        job = Mock(Job)
        job.metabolites2csv.return_value = csv
        request = testing.DummyRequest(params={
                                               'start':0,
                                               'limit':10
                                               })
        views = Views(request)
        views.job = Mock(return_value=job)
        views.metabolitesjson = Mock(return_value={
                                                   'rows':[]
                                                   })

        response = views.metabolitescsv()

        self.assertEqual(response.content_type, 'text/csv')
        self.assertEqual(response.body, 'bla')

    def test_metabolitessdf(self):
        job = Mock(Job)
        job.metabolites2sdf.return_value = 'bla'
        request = testing.DummyRequest(params={
                                               'start':0,
                                               'limit':10
                                               })
        views = Views(request)
        views.job = Mock(return_value=job)
        views.metabolitesjson = Mock(return_value={
                                                   'rows':[]
                                                   })
        response = views.metabolitessdf()

        self.assertEqual(response.content_type, 'chemical/x-mdl-sdfile')
        self.assertEqual(response.body, 'bla')


    def test_chromatogramjson(self):
        request = testing.DummyRequest()
        views = Views(request)
        views.job = Mock(return_value=self.fake_job())

        response  = views.chromatogramjson()

        self.assertEqual(response, [1,2,3])

    def test_mspectrajson(self):
        request = testing.DummyRequest(matchdict={'scanid':641})
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        views.mspectrajson()

        job.mspectra.assert_called_with(641, None)

    def test_mspectrajson_withmslevel(self):
        request = testing.DummyRequest(
                                       matchdict={'scanid':641},
                                       params={'mslevel':3}
                                       )
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        views.mspectrajson()

        job.mspectra.assert_called_with(641, 3)

    def test_mspectra_withoutscanid_notfound(self):
        request = testing.DummyRequest(matchdict={'scanid':641})
        views = Views(request)
        job = self.fake_job()
        from magmaweb.job import ScanNotFound
        job.mspectra.side_effect = ScanNotFound()
        views.job = Mock(return_value=job)

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            views.mspectrajson()

    def test_extractedionchromatogram(self):
        request = testing.DummyRequest(matchdict={'metid':72})
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        response = views.extractedionchromatogram()

        self.assertEqual(response, {
            'chromatogram': [1,2,3], 'scans': [4,5]
        })
        job.extractedIonChromatogram.assert_called_with(72)
        job.scansWithMetabolites.assert_called_with(metid=72)


    def test_fragments_metabolitewithoutfragments(self):
        request = testing.DummyRequest(matchdict={
                                                  'metid':72,
                                                  'scanid':641
                                                  },
                                       params={'node':''})
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        views.fragments()

        job.fragments.assert_called_with(metid=72, scanid=641, node='')

    def test_fragments_filteronmslevel(self):
        request = testing.DummyRequest(matchdict={
                                                  'metid':72,
                                                  'scanid':641
                                                  },
                                       params={'node':1709})
        views = Views(request)
        job = self.fake_job()
        views.job = Mock(return_value=job)

        views.fragments()

        job.fragments.assert_called_with(metid=72, scanid=641, node=1709)

    def test_fragments_badmetabolite_notfound(self):
        request = testing.DummyRequest(matchdict={
                                                  'metid':72,
                                                  'scanid':641
                                                  },
                                       params={'node':''})
        views = Views(request)
        job = self.fake_job()
        from magmaweb.job import FragmentNotFound
        job.fragments.side_effect = FragmentNotFound
        views.job = Mock(return_value=job)

        from pyramid.httpexceptions import HTTPNotFound
        with self.assertRaises(HTTPNotFound):
            views.fragments()

    def test_stderr(self):
        request = testing.DummyRequest()
        views = Views(request)
        job = self.fake_job()
        import StringIO
        log = StringIO.StringIO()
        log.write('bla')
        log.seek(0)
        job.stderr.return_value = log
        views.job = Mock(return_value=job)

        response = views.stderr()

        self.assertEqual(response.content_type, 'text/plain')
        self.assertMultiLineEqual(response.app_iter.read(), 'bla')

    def test_runinfojson(self):
        request = testing.DummyRequest()
        views = Views(request)
        job = self.fake_job()
        from magmaweb.models import Run
        job.runInfo.return_value = Run(
            n_reaction_steps=2, metabolism_types='phase1,phase2' ,
            ionisation_mode=-1, skip_fragmentation=True,
            ms_intensity_cutoff=200000.0, msms_intensity_cutoff=0.5,
            mz_precision=0.01, use_all_peaks=True,
            ms_filename = 'F123456.mzxml', abs_peak_cutoff=1000,
            rel_peak_cutoff=0.001, max_ms_level=3, precursor_mz_precision=0.01,
            max_broken_bonds=4, description='My first description'
        )
        views.job = Mock(return_value=job)

        response = views.runinfojson()

        self.assertEqual(response, {
                                    'success': True,
                                    'data': dict(
                                                 n_reaction_steps=2,
                                                 metabolism_types=['phase1', 'phase2'],
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
        self.maxDiff =200000
        request = testing.DummyRequest()
        views = Views(request)
        job = self.fake_job()
        from magmaweb.models import Run
        job.runInfo.return_value = Run(
            abs_peak_cutoff=1100,
            rel_peak_cutoff=0.012
        )
        views.job = Mock(return_value=job)

        response = views.runinfojson()

        self.assertEqual(response, {
                                    'success': True,
                                    'data': dict(
                                                 n_reaction_steps=2,
                                                 metabolism_types=['phase1', 'phase2'],
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
        views = Views(request)
        job = self.fake_job()
        job.runInfo.return_value = None
        views.job = Mock(return_value=job)

        response = views.runinfojson()

        self.assertEqual(response, {
                                    'success': True,
                                    'data': dict(
                                                 n_reaction_steps=2,
                                                 metabolism_types=['phase1', 'phase2'],
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


    def test_defaultsjson(self):
        request = testing.DummyRequest()
        views = Views(request)
        response = views.defaults()

        self.assertEqual(response, {
                                    'success': True,
                                    'data': dict(
                                                 n_reaction_steps=2,
                                                 metabolism_types=['phase1', 'phase2'],
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

