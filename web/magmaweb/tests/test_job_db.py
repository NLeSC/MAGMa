"""Tests for magmaweb.job.JobDb"""
import unittest
from magmaweb.job import JobDb
from magmaweb.models import Scan, Peak, Run, Metabolite, Fragment
from magmaweb.tests.test_job import initTestingDB


class JobDbTestCaseAbstract(unittest.TestCase):
    def setUp(self):
        # mock job session
        self.session = initTestingDB()
        self.job = JobDb(self.session)


class JobDbTestCase(JobDbTestCaseAbstract):
    def test_construct(self):
        self.assertEqual(self.job.session, self.session)

    def test_runInfo(self):
        runInfo = self.job.runInfo()
        self.assertEqual(runInfo.n_reaction_steps, 2)
        self.assertEqual(runInfo.metabolism_types, u"phase1,phase2")

        self.assertEqual(runInfo.ms_filename, u'F123456.mzxml')
        self.assertEqual(runInfo.abs_peak_cutoff, 1000)
        self.assertEqual(runInfo.max_ms_level, 3)
        self.assertEqual(runInfo.precursor_mz_precision, 10)

        self.assertEqual(runInfo.ionisation_mode, -1)
        self.assertEqual(runInfo.max_broken_bonds, 4)
        self.assertEqual(runInfo.ms_intensity_cutoff, 200000.0)
        self.assertEqual(runInfo.msms_intensity_cutoff, 50)
        self.assertEqual(runInfo.mz_precision, 10)
        self.assertEqual(runInfo.mz_precision_abs, 0.002)
        self.assertEqual(runInfo.description, u'My first description')
        self.assertEqual(runInfo.max_water_losses, 1)

    def test_runInfo_maxrunid(self):
        self.session.add(Run(
            n_reaction_steps=2, metabolism_types=u'phase1,phase2',
            ionisation_mode=-1, skip_fragmentation=True,
            ms_intensity_cutoff=200000.0, msms_intensity_cutoff=50,
            mz_precision=10, mz_precision_abs=0.002, use_all_peaks=True,
            ms_filename=u'F123456.mzxml', abs_peak_cutoff=1000,
            max_ms_level=3, precursor_mz_precision=10,
            max_broken_bonds=4, description=u'My second description',
            max_water_losses=1,
        ))

        runInfo = self.job.runInfo()

        # run with highest id is returned
        self.assertEqual(runInfo.description, u'My second description')

    def test_maxMSLevel(self):
        maxmslevel = self.job.maxMSLevel()
        self.assertEqual(maxmslevel, 3)

    def test_maxMSLevel_withoutscans(self):
        self.session.query(Peak).delete()
        self.session.query(Scan).delete()
        maxmslevel = self.job.maxMSLevel()
        self.assertEqual(maxmslevel, 0)

    def test_extractedionchromatogram(self):
        metid = 72
        eic = self.job.extractedIonChromatogram(metid)
        self.assertEqual(eic, [{
            'rt': 933.317,
            'intensity': 345608.65625
        }, {
            'rt': 1254.15,
            'intensity': 0
        }])

    def test_chromatogram(self):
        expected_chromatogram = {'scans': [{'id': 641, 'rt': 933.317,
                                            'intensity': 807577.0, 'ap': 0},
                                           {'id': 870, 'rt': 1254.15,
                                            'intensity': 1972180.0, 'ap': 0},
                                           ],
                                 'cutoff': 200000.0,
                                 }
        self.assertEqual(self.job.chromatogram(), expected_chromatogram)

    def test_chromatogram_with_assigned_peaks(self):
        metid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz, metid)

        expected_chromatogram = {'scans': [{'id': 641, 'rt': 933.317,
                                            'intensity': 807577.0, 'ap': 1},
                                           {'id': 870, 'rt': 1254.15,
                                            'intensity': 1972180.0, 'ap': 0},
                                           ],
                                 'cutoff': 200000.0,
                                 }
        self.assertEqual(self.job.chromatogram(), expected_chromatogram)

    def test_metabolitesTotalCount(self):
        self.assertEqual(self.job.metabolitesTotalCount(), 2)

    def test_assign_metabolite2peak(self):
        metid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz, metid)

        q = self.session.query(Peak.assigned_metid)
        q = q.filter(Peak.scanid == scanid)
        expected_metid = q.filter(Peak.mz == mz).scalar()
        self.assertEqual(expected_metid, metid)

    def test_assign_metabolite2peak_withoffset(self):
        metid = 72
        scanid = 641
        offset = 8e-7
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz + offset, metid)

        q = self.session.query(Peak.assigned_metid)
        q = q.filter(Peak.scanid == scanid)
        expected_metid = q.filter(Peak.mz == mz).scalar()
        self.assertEqual(expected_metid, metid)

    def test_unassign_metabolite2peak(self):
        metid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz, metid)

        self.job.unassign_metabolite2peak(scanid, mz)

        q = self.session.query(Peak.assigned_metid)
        q = q.filter(Peak.scanid == scanid)
        self.assertIsNone(q.filter(Peak.mz == mz).scalar())

    def test_unassign_metabolite2peak_withoffset(self):
        metid = 72
        scanid = 641
        offset = 8e-7
        mz = 109.0295639038086
        self.job.assign_metabolite2peak(scanid, mz, metid)

        self.job.unassign_metabolite2peak(scanid, mz + offset)

        q = self.session.query(Peak.assigned_metid)
        q = q.filter(Peak.scanid == scanid)
        self.assertIsNone(q.filter(Peak.mz == mz).scalar())

    def test_hasMolecules(self):
        self.assertTrue(self.job.hasMolecules())

    def test_hasNoMolecules(self):
        self.job.session.query(Metabolite).delete()
        self.assertFalse(self.job.hasMolecules())

    def test_hasMspectras(self):
        self.assertTrue(self.job.hasMspectras())

    def test_hasNoMspectras(self):
        self.job.session.query(Peak).delete()
        self.assertFalse(self.job.hasMspectras())

    def test_hasFragments(self):
        self.assertTrue(self.job.hasFragments())

    def test_hasNoFragments(self):
        self.job.session.query(Fragment).delete()
        self.assertFalse(self.job.hasFragments())


class JobDbEmptyDatasetTestCase(unittest.TestCase):
    def setUp(self):
        # mock job session
        self.session = initTestingDB(dataset=None)
        self.job = JobDb(self.session)

    def test_runInfo(self):
        runInfo = self.job.runInfo()
        self.assertIsNone(runInfo)

    def test_maxMSLevel(self):
        maxmslevel = self.job.maxMSLevel()
        self.assertEqual(maxmslevel, 0)

    def test_chromatogram(self):
        expected_chromatogram = {'scans': [], 'cutoff': None}
        self.assertEqual(self.job.chromatogram(), expected_chromatogram)

    def test_metabolitesTotalCount(self):
        self.assertEqual(self.job.metabolitesTotalCount(), 0)


class JobDbMetabolitesTestCase(JobDbTestCaseAbstract):
    def test_default(self):
        response = self.job.metabolites()
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += '?cid=152432">CID: 152432</a>'
        self.assertEquals(
            response,
            {
                'total': 2,
                'rows': [{
                    'metid': 72,
                    'isquery': True,
                    'level': 0,
                    'mol': u'Molfile',
                    'molformula': u'C6H6O2',
                    'nhits': 1,
                    'origin': u'pyrocatechol',
                    'probability': 1.0,
                    'reactionsequence': [u'PARENT'],
                    'smiles': u'Oc1ccccc1O',
                    'mim': 110.03677, 'logp':1.231,
                    'assigned': False,
                    'reference': url1
                }, {
                    'isquery': True, 'level': 0, 'metid': 352,
                    'mol': u"Molfile of dihydroxyphenyl-valerolactone",
                    'molformula': u"C11H12O4",
                    'nhits': 1,
                    'origin': u"dihydroxyphenyl-valerolactone",
                    'probability': 1,
                    'reactionsequence': [u"PARENT", u"CHILD"],
                    'smiles': u"O=C1OC(Cc2ccc(O)c(O)c2)CC1",
                    'mim': 208.07355, 'logp':2.763,
                    'assigned': False,
                    'reference': url2
                }]
            }
        )

    def test_scanid(self):
        response = self.job.metabolites(scanid=641)
        self.assertIn('score', response['rows'][0])
        self.assertIn('deltappm', response['rows'][0])
        self.assertEqual(response['total'], 1)

    def test_filteredon_nrscanseq(self):
        response = self.job.metabolites(filters=[{"type": "numeric",
                                                  "comparison": "eq",
                                                  "value": 1,
                                                  "field": "nhits"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscansgt(self):
        response = self.job.metabolites(filters=[{"type": "numeric",
                                                  "comparison": "lt",
                                                  "value": 2,
                                                  "field": "nhits"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscanslt(self):
        response = self.job.metabolites(filters=[{"type": "numeric",
                                                  "comparison": "gt",
                                                  "value": 0,
                                                  "field": "nhits"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_isquery(self):
        response = self.job.metabolites(filters=[{"type": "boolean",
                                                  "value": True,
                                                  "field": "isquery"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_molformula(self):
        response = self.job.metabolites(filters=[{"type": "string",
                                                  "value": u"C6",
                                                  "field": "molformula"}])

        self.assertEqual(response['total'], 1)

    def test_filteredon_level(self):
        response = self.job.metabolites(filters=[{"type": "list",
                                                  "value": [0, 1, 2],
                                                  "field": "level"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_score(self):
        filters = [{"type": "numeric", "comparison": "eq",
                    "value": 200, "field": "score"}]

        response = self.job.metabolites(scanid=641, filters=filters)

        self.assertEqual(response['total'], 1)

    def test_filteredon_score_without_scan(self):
        from magmaweb.job import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(filters=[{"type": "numeric",
                                           "comparison": "eq",
                                           "value": 200,
                                           "field": "score"}])

    def test_filteredon_deltappm(self):
        filters = [{"type": "numeric", "comparison": "eq",
                    "value": -1.84815979523607e-08, "field": "deltappm"}]

        response = self.job.metabolites(scanid=641, filters=filters)

        self.assertEqual(response['total'], 1)

    def test_filteredon_deltappm_without_scan(self):
        from magmaweb.job import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(filters=[{"type": "numeric",
                                           "comparison": "eq",
                                           "value": -1.84815979523607e-08,
                                           "field": "deltappm"}])

    def test_filteredon_not_assigned(self):
        response = self.job.metabolites(filters=[{"type": "boolean",
                                                  "value": False,
                                                  "field": "assigned"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_assigned(self):
        response = self.job.metabolites(filters=[{"type": "boolean",
                                                  "value": True,
                                                  "field": "assigned"}])

        self.assertEqual(response['total'], 0)

    def test_sort_probmet(self):
        response = self.job.metabolites(sorts=[{"property": "probability",
                                                "direction": "DESC"},
                                               {"property": "metid",
                                                "direction": "ASC"}])

        self.assertEqual(response['total'], 2)

    def test_sort_nrscans(self):
        response = self.job.metabolites(sorts=[{"property": "nhits",
                                                "direction": "DESC"}])

        self.assertEqual(response['total'], 2)

    def test_sort_assigned(self):
        response = self.job.metabolites(sorts=[{"property": "assigned",
                                                "direction": "DESC"}])

        self.assertEqual(response['total'], 2)

    def test_sort_score(self):
        sorts = [{"property": "score", "direction": "DESC"}]

        response = self.job.metabolites(scanid=641, sorts=sorts)

        self.assertEqual(response['total'], 1)

    def test_sort_score_without_scan(self):
        from magmaweb.job import ScanRequiredError
        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(sorts=[{"property": "score",
                                         "direction": "DESC"}])

    def test_sort_deltappm(self):
        sorts = [{"property": "deltappm", "direction": "DESC"}]

        response = self.job.metabolites(scanid=641, sorts=sorts)

        self.assertEqual(response['total'], 1)

    def test_sort_deltappm_without_scan(self):
        sorts = [{"property": "deltappm", "direction": "DESC"}]
        from magmaweb.job import ScanRequiredError

        with self.assertRaises(ScanRequiredError):
            self.job.metabolites(sorts=sorts)


class JobDbMetabolites2csvTestCase(JobDbTestCaseAbstract):
    def test_it(self):
        csvfile = self.job.metabolites2csv(self.job.metabolites()['rows'])
        import csv
        import StringIO
        expected_csvfile = StringIO.StringIO()
        cols = ['origin', 'smiles', 'probability', 'reactionsequence',
                'nhits', 'molformula', 'mim', 'isquery', 'logp', 'reference']
        csvwriter = csv.DictWriter(expected_csvfile, cols)
        csvwriter.writeheader()
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += '?cid=152432">CID: 152432</a>'
        csvwriter.writerow({'origin': 'pyrocatechol',
                            'smiles': 'Oc1ccccc1O',
                            'probability': 1.0,
                            'reactionsequence': 'PARENT',
                            'nhits': 1,
                            'molformula': 'C6H6O2',
                            'isquery': True,
                            'mim': 110.03677,
                            'logp': 1.231,
                            'reference': url1,
                            })
        csvwriter.writerow({'origin': 'dihydroxyphenyl-valerolactone',
                            'smiles': 'O=C1OC(Cc2ccc(O)c(O)c2)CC1',
                            'probability': 1.0,
                            'reactionsequence': 'PARENT|CHILD',
                            'nhits': 1,
                            'molformula': 'C11H12O4',
                            'isquery': True,
                            'mim': 208.07355,
                            'logp': 2.763,
                            'reference': url2,
                            })
        self.assertMultiLineEqual(csvfile.getvalue(),
                                  expected_csvfile.getvalue())

    def test_with_score(self):
        mets = self.job.metabolites(scanid=641)['rows']
        csvfile = self.job.metabolites2csv(mets)
        import csv
        import StringIO
        expected_csvfile = StringIO.StringIO()
        cols = ['origin', 'smiles', 'probability', 'reactionsequence',
                'nhits', 'molformula', 'mim', 'isquery', 'logp', 'reference',
                'score']
        csvwriter = csv.DictWriter(expected_csvfile, cols)
        csvwriter.writeheader()
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        csvwriter.writerow({'origin': 'pyrocatechol', 'smiles': 'Oc1ccccc1O',
                            'probability': 1.0, 'reactionsequence': 'PARENT',
                            'nhits': 1, 'molformula': 'C6H6O2',
                            'isquery': True, 'score': 200.0, 'mim': 110.03677,
                            'logp': 1.231,
                            'reference': url1
                            })
        self.assertMultiLineEqual(csvfile.getvalue(),
                                  expected_csvfile.getvalue())

    def test_some_columns(self):
        cols = ['origin', 'mim']
        mets = self.job.metabolites(scanid=641)['rows']

        response = self.job.metabolites2csv(mets, cols=cols)

        self.assertEquals(
            response.getvalue(),
            'origin,mim\r\n' +
            'pyrocatechol,110.03677\r\n'
        )


class JobMetabolites2sdfTestCase(JobDbTestCaseAbstract):
    def test_it(self):
        sdffile = self.job.metabolites2sdf(self.job.metabolites()['rows'])
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += '?cid=152432">CID: 152432</a>'

        expected_sdf = """Molfile> <origin>
pyrocatechol

> <smiles>
Oc1ccccc1O

> <probability>
1.0

> <reactionsequence>
PARENT

> <nhits>
1

> <molformula>
C6H6O2

> <mim>
110.03677

> <logp>
1.231

> <reference>
{url1}

$$$$
Molfile of dihydroxyphenyl-valerolactone> <origin>
dihydroxyphenyl-valerolactone

> <smiles>
O=C1OC(Cc2ccc(O)c(O)c2)CC1

> <probability>
1.0

> <reactionsequence>
PARENT
CHILD

> <nhits>
1

> <molformula>
C11H12O4

> <mim>
208.07355

> <logp>
2.763

> <reference>
{url2}

$$$$
""".format(url1=url1, url2=url2)

        self.assertMultiLineEqual(sdffile, expected_sdf)

    def test_with_scan(self):
        """Include score prop"""
        mets = self.job.metabolites(scanid=641)['rows']
        sdffile = self.job.metabolites2sdf(mets)
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        expected_sdf = """Molfile> <origin>
pyrocatechol

> <smiles>
Oc1ccccc1O

> <probability>
1.0

> <reactionsequence>
PARENT

> <nhits>
1

> <molformula>
C6H6O2

> <mim>
110.03677

> <logp>
1.231

> <reference>
{url1}

> <score>
200.0

$$$$
""".format(url1=url1)

        self.assertMultiLineEqual(sdffile, expected_sdf)

    def test_some_columns(self):
        cols = ['origin', 'mim']

        mets = self.job.metabolites(scanid=641)['rows']
        sdffile = self.job.metabolites2sdf(mets, cols=cols)

        expected_sdf = """Molfile> <origin>
pyrocatechol

> <mim>
110.03677

$$$$
"""

        self.assertMultiLineEqual(sdffile, expected_sdf)


class JobScansWithMetabolitesTestCase(JobDbTestCaseAbstract):
    def test_metid(self):
        response = self.job.scansWithMetabolites(metid=72)
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_all(self):
        response = self.job.scansWithMetabolites()
        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_molformula(self):
        filters = [{"type": "string", "value": u"C6", "field": "molformula"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_nrscans(self):
        filters = [{"type": "numeric", "value": "1",
                    "comparison": "eq", "field": "nhits"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_score(self):
        filters = [{"type": "numeric", "value": "200",
                    "comparison": "eq", "field": "score"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317}
        ])

    def test_filteredon_not_assigned(self):
        filters = [{"type": "boolean", "value": False, "field": "assigned"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_filteredon_assigned(self):
        filters = [{"type": "boolean", "value": True, "field": "assigned"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [])

    def test_filteredon_nrscansgt_and_molformula(self):
        filters = [{"type": "numeric", "comparison": "lt",
                    "value": 2, "field": "nhits"},
                   {"type": "string", "value": u"C6", "field": "molformula"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(len(response), 1)

    def test_filteron_deltappm(self):
        filters = [{"type": "numeric", "comparison": "gt",
                    "value": 0, "field": "nhits"},
                   {"type": "numeric", "value": "200",
                    "comparison": "lt", "field": "deltappm"}]

        response = self.job.scansWithMetabolites(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])


class JobMSpectraTestCase(JobDbTestCaseAbstract):
    def test_scanonly(self):
        self.assertEqual(
            self.job.mspectra(641),
            {
                'peaks': [
                    {'intensity': 345608.65625,
                     'mz': 109.0295639038086,
                     'assigned_metid': None},
                    {'intensity': 807576.625,
                     'mz': 305.033508300781,
                     'assigned_metid': None}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': {'id': None, 'mz': None}
            }
        )

    def test_withmslevel(self):
        self.assertEqual(
            self.job.mspectra(641, 1),
            {
                'peaks': [
                    {'intensity': 345608.65625,
                     'mz': 109.0295639038086,
                     'assigned_metid': None},
                    {'intensity': 807576.625,
                     'mz': 305.033508300781,
                     'assigned_metid': None}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': {'id': None, 'mz': None}
            }
        )

    def test_notfound(self):
        from magmaweb.job import ScanNotFound
        with self.assertRaises(ScanNotFound):
            self.job.mspectra(123)

    def test_lvl2scan(self):
        response = self.job.mspectra(871, 2)
        self.assertEqual(response, {
            'peaks': [
                {'intensity': 211603.046875,
                 'mz': 123.04508972168,
                 'assigned_metid': None},
                {'intensity': 279010.28125,
                 'mz': 163.076232910156,
                 'assigned_metid': None}
            ],
            'cutoff': 13950500.0,
            'mslevel': 2,
            'precursor': {'id': 870, 'mz': 207.0663147}
        })

    def test_withassigned_met2peak(self):
        self.job.assign_metabolite2peak(641, 109.0295639038086, 72)
        self.assertEqual(
            self.job.mspectra(641),
            {
                'peaks': [
                    {'intensity': 345608.65625,
                     'mz': 109.0295639038086,
                     'assigned_metid': 72},
                    {'intensity': 807576.625,
                     'mz': 305.033508300781,
                     'assigned_metid': None}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': {'id': None, 'mz': None}
            }
        )


class JobFragmentsTestCase(JobDbTestCaseAbstract):
    def test_metabolitewithoutfragments(self):
        self.maxDiff = None
        response = self.job.fragments(metid=72, scanid=641, node='root')
        self.assertEqual(response, {
            'children': [{
                'atoms': u'0,1,2,3,4,5,6,7',
                'children': [],
                'deltah': -1.0,
                'deltappm': -1.84815979523607e-08,
                'expanded': True,
                'fragid': 948,
                'leaf': True,
                'mass': 110.0367794368,
                'metid': 72,
                'mol': u'Molfile',
                'formula': u'C5H4',
                'mslevel': 1,
                'mz': 109.0295639038086,
                'scanid': 641,
                'score': 200.0,
                'isAssigned': False,
            }], 'expanded': True
        })

    def test_metabolitewithfragments(self):
        self.maxDiff = None
        response = self.job.fragments(metid=352, scanid=870, node='root')
        self.assertEqual(response, {
            'children': [{
                'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14',
                'children': [{
                    'atoms': "6,7,8,9,10,11,12,13,14",
                    'deltah': 0,
                    'deltappm': 3.943698856902144e-12,
                    'expanded': True,
                    'fragid': 1708,
                    'leaf': True,
                    'mass': 123.0446044689,
                    'metid': 352,
                    'mol': "Molfile of dihydroxyphenyl-valerolactone",
                    'mslevel': 2,
                    'mz': 123.04508972167969,
                    'scanid': 871,
                    'score': 201,
                    'formula': 'C3H5O3',
                }, {
                    'atoms': "3,4,5,6,7,8,9,10,11,12,13,14",
                    'deltah': -1,
                    'deltappm': -1.235815738001507e-08,
                    'expanded': False,
                    'fragid': 1709,
                    'leaf': False,
                    'mass': 164.08372962939995,
                    'metid': 352,
                    'mol': "Molfile of dihydroxyphenyl-valerolactone",
                    'mslevel': 2,
                    'mz': 163.07623291015625,
                    'scanid': 871,
                    'score': 65,
                    'formula': 'C3H5O3',
                }],
                'deltah': -1,
                'deltappm': -8.18675317722029e-09,
                'expanded': True,
                'fragid': 1707,
                'leaf': False,
                'mass': 208.0735588736,
                'metid': 352,
                'mol': u'Molfile of dihydroxyphenyl-valerolactone',
                'mslevel': 1,
                'mz': 207.066284179688,
                'scanid': 870,
                'score': 100,
                'isAssigned': False,
                'formula': 'C3H5O3',
            }], 'expanded': True
        })

    def test_lvl3fragments(self):
        response = self.job.fragments(metid=352, scanid=870, node=1709)
        self.assertEqual(response, [{
            'atoms': "4,5,6,7,8,9,11,13,14",
            'deltah': 3,
            'expanded': True,
            'fragid': 1710,
            'leaf': True,
            'mass': 116.0626002568,
            'metid': 352,
            'mol': "Molfile of dihydroxyphenyl-valerolactone",
            'mslevel': 3,
            'mz': 119.08654022216797,
            'scanid': 872,
            'score': 4,
            'deltappm': 5.0781684060061766e-08,
            'formula': 'C3H5O3',
        }])

    def test_metabolitewithassignedpeak(self):
        self.job.assign_metabolite2peak(641, 109.0295639038086, 72)
        response = self.job.fragments(metid=72, scanid=641, node='root')
        self.assertEqual(response, {
            'children': [{
                'atoms': u'0,1,2,3,4,5,6,7',
                'children': [],
                'deltah': -1.0,
                'deltappm': -1.84815979523607e-08,
                'expanded': True,
                'fragid': 948,
                'leaf': True,
                'mass': 110.0367794368,
                'metid': 72,
                'mol': u'Molfile',
                'mslevel': 1,
                'mz': 109.0295639038086,
                'scanid': 641,
                'score': 200.0,
                'isAssigned': True,
                'formula': "C5H4",
            }], 'expanded': True
        })

    def test_badfragment(self):
        from magmaweb.job import FragmentNotFound
        with self.assertRaises(FragmentNotFound):
            self.job.fragments(metid=70002, scanid=641, node='root')


class JobWithAllPeaksTestCase(unittest.TestCase):
    def setUp(self):
        self.job = JobDb(initTestingDB(dataset='useallpeaks'))

    def test_default(self):
        response = self.job.metabolites()
        self.assertEquals(
            response,
            {
                'total': 1,
                'rows': [{
                    'metid': 12,
                    'isquery': False,
                    'level': 1,
                    'mol': u'Molfile',
                    'molformula': u'C11H12O7S',
                    'nhits': 1,
                    'origin': u'5-(3,4)-dihydroxyphenyl-g-valerolactone (F)',
                    'probability': 0.119004,
                    'reactionsequence': [u'sulfation_(aromatic_hydroxyl)'],
                    'smiles': u'Oc1ccc(CC2OC(=O)CC2)cc1OS(O)(=O)=O',
                    'mim': 288.0303734299, 'logp':1.9027,
                    'assigned': False,
                    'reference': ''
                }]
            }
        )

    def test_lvl1fragments(self):
        response = self.job.fragments(metid=12, scanid=1, node='root')
        self.assertEqual(response, {
            'children': [{
                'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
                'children': [],
                'deltah': -1.0,
                'deltappm': -7.046696487857745e-09,
                'expanded': True,
                'fragid': 17,
                'leaf': True,
                'mass': 288.0303734299,
                'metid': 12,
                'mol': u'Molfile',
                'mslevel': 1,
                'mz': 287.015686035156,
                'scanid': 1,
                'score': 0.0,
                'isAssigned': False,
                'formula': "C4H6O2",
            }, {
                'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
                'deltah': -1.0,
                'deltappm': -7.020570507205176e-09,
                'expanded': True,
                'fragid': 18,
                'leaf': False,
                'mz': 287.023132324219,
                'metid': 12,
                'mol': u'Molfile',
                'mslevel': 1,
                'mass': 288.0303734299,
                'scanid': 1,
                'score': 0.5,
                'isAssigned': False,
                'formula': "C4H6O2",
                'children': [{
                    'fragid': 19,
                    'metid': 12,
                    'scanid': 2,
                    'mz': 207.066223144531,
                    'mass': 207.0657338415,
                    'score': 0.5,
                    'mol': u'Molfile',
                    'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14',
                    'expanded': True,
                    'leaf': True,
                    'mslevel': 2,
                    'deltah': 0.0,
                    'deltappm': 2.3630267823129625e-12,
                    'formula': "C4H6O2",
                }, {
                    'fragid': 20,
                    'metid': 12,
                    'scanid': 2,
                    'mz': 287.022827148438,
                    'mass': 288.0303734299,
                    'score': 0.0,
                    'mol': u'Molfile',
                    'atoms': u'0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18',
                    'expanded': True,
                    'leaf': True,
                    'mslevel': 2,
                    'deltah': -1.0,
                    'deltappm': -7.021641217476223e-09,
                    'formula': "C4H6O2",
                }]
            }], 'expanded': True
        })
