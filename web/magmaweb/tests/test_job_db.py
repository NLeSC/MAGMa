"""Tests for magmaweb.job.JobDb"""
import unittest
from magmaweb.job import JobDb, ScanRequiredError
from magmaweb.models import Scan, Peak, Run, Molecule, Fragment
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
        molid = 72
        eic = self.job.extractedIonChromatogram(molid)
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
        molid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_molecule2peak(scanid, mz, molid)

        expected_chromatogram = {'scans': [{'id': 641, 'rt': 933.317,
                                            'intensity': 807577.0, 'ap': 1},
                                           {'id': 870, 'rt': 1254.15,
                                            'intensity': 1972180.0, 'ap': 0},
                                           ],
                                 'cutoff': 200000.0,
                                 }
        self.assertEqual(self.job.chromatogram(), expected_chromatogram)

    def test_moleculesTotalCount(self):
        self.assertEqual(self.job.moleculesTotalCount(), 2)

    def test_assign_molecule2peak(self):
        molid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_molecule2peak(scanid, mz, molid)

        q = self.session.query(Peak.assigned_molid)
        q = q.filter(Peak.scanid == scanid)
        expected_molid = q.filter(Peak.mz == mz).scalar()
        self.assertEqual(expected_molid, molid)

    def test_assign_molecule2peak_withoffset(self):
        molid = 72
        scanid = 641
        offset = 8e-7
        mz = 109.0295639038086
        self.job.assign_molecule2peak(scanid, mz + offset, molid)

        q = self.session.query(Peak.assigned_molid)
        q = q.filter(Peak.scanid == scanid)
        expected_molid = q.filter(Peak.mz == mz).scalar()
        self.assertEqual(expected_molid, molid)

    def test_unassign_molecule2peak(self):
        molid = 72
        scanid = 641
        mz = 109.0295639038086
        self.job.assign_molecule2peak(scanid, mz, molid)

        self.job.unassign_molecule2peak(scanid, mz)

        q = self.session.query(Peak.assigned_molid)
        q = q.filter(Peak.scanid == scanid)
        self.assertIsNone(q.filter(Peak.mz == mz).scalar())

    def test_unassign_molecule2peak_withoffset(self):
        molid = 72
        scanid = 641
        offset = 8e-7
        mz = 109.0295639038086
        self.job.assign_molecule2peak(scanid, mz, molid)

        self.job.unassign_molecule2peak(scanid, mz + offset)

        q = self.session.query(Peak.assigned_molid)
        q = q.filter(Peak.scanid == scanid)
        self.assertIsNone(q.filter(Peak.mz == mz).scalar())

    def test_hasMolecules(self):
        self.assertTrue(self.job.hasMolecules())

    def test_hasNoMolecules(self):
        self.job.session.query(Molecule).delete()
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

    def test_moleculesTotalCount(self):
        self.assertEqual(self.job.moleculesTotalCount(), 0)


class JobDbMoleculesTestCase(JobDbTestCaseAbstract):

    def test_default(self):
        response = self.job.molecules()
        url1 = u'<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += u'?cid=289">CID: 289</a>'
        url2 = u'<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += u'?cid=152432">CID: 152432</a>'
        self.assertEquals(
            response,
            {
                'total': 2,
                'page': 1,
                'rows': [{
                    'molid': 72,
                    'predicted': False,
                    'mol': u'Molfile',
                    'formula': u'C6H6O2',
                    'nhits': 1,
                    'name': u'pyrocatechol',
                    'refscore': 1.0,
                    'reactionsequence': {
                        u'reactantof': {
                            u'esterase': {
                                u'nr': 2,
                                u'nrp': 1
                            }
                        }
                    },
                    'smiles': u'C1=CC=C(C(=C1)O)O',
                    'inchikey14': u'YCIMNLLNPGFGHC',
                    'mim': 110.03677, 'logp': 1.231,
                    'assigned': False,
                    'reference': url1
                }, {
                    'predicted': False, 'molid': 352,
                    'mol': u"Molfile of dihydroxyphenyl-valerolactone",
                    'formula': u"C11H12O4",
                    'nhits': 1,
                    'name': u"dihydroxyphenyl-valerolactone",
                    'refscore': 1.0,
                    'reactionsequence': {
                        u'productof': {
                            u'theogallin': {
                                u'nr': 1,
                                u'nrp': 0
                            }
                        }
                    },
                    'smiles': u"O=C1CCC(Cc2ccc(O)c(O)c2)O1",
                    'inchikey14': u'ZNXXWTPQHVLMQT',
                    'mim': 208.07355, 'logp': 2.763,
                    'assigned': False,
                    'reference': url2
                }]
            }
        )

    def test_scanid(self):
        response = self.job.molecules(scanid=641)
        self.assertIn('score', response['rows'][0])
        self.assertIn('deltappm', response['rows'][0])
        self.assertIn('mz', response['rows'][0])
        self.assertEqual(response['total'], 1)

    def test_filteredon_nrscanseq(self):
        response = self.job.molecules(filters=[{"type": "numeric",
                                                "comparison": "eq",
                                                "value": 1,
                                                "field": "nhits"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscansgt(self):
        response = self.job.molecules(filters=[{"type": "numeric",
                                                "comparison": "lt",
                                                "value": 2,
                                                "field": "nhits"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_nrscanslt(self):
        response = self.job.molecules(filters=[{"type": "numeric",
                                                "comparison": "gt",
                                                "value": 0,
                                                "field": "nhits"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_predicted(self):
        response = self.job.molecules(filters=[{"type": "boolean",
                                                "value": False,
                                                "field": "predicted"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_formula(self):
        response = self.job.molecules(filters=[{"type": "string",
                                                "value": u"C6",
                                                "field": "formula"}])

        self.assertEqual(response['total'], 1)

    def test_filteredon_score(self):
        filters = [{"type": "numeric", "comparison": "eq",
                    "value": 200, "field": "score"}]

        response = self.job.molecules(scanid=641, filters=filters)

        self.assertEqual(response['total'], 1)

    def test_filteredon_score_without_scan(self):
        with self.assertRaises(ScanRequiredError):
            self.job.molecules(filters=[{"type": "numeric",
                                         "comparison": "eq",
                                         "value": 200,
                                         "field": "score"}])

    def test_filteredon_deltappm(self):
        filters = [{"type": "numeric", "comparison": "eq",
                    "value": -1.84815979523607e-08, "field": "deltappm"}]

        response = self.job.molecules(scanid=641, filters=filters)

        self.assertEqual(response['total'], 1)

    def test_filteredon_deltappm_without_scan(self):
        with self.assertRaises(ScanRequiredError):
            self.job.molecules(filters=[{"type": "numeric",
                                         "comparison": "eq",
                                         "value": -1.84815979523607e-08,
                                         "field": "deltappm"}])

    def test_filteredon_mz(self):
        response = self.job.molecules(scanid=641, mz=109.0295639038086)

        self.assertEqual(response['total'], 1)

    def test_filteredon_mz_without_scan(self):
        with self.assertRaises(ScanRequiredError):
            self.job.molecules(mz=109.0295639038086)

    def test_filteredon_not_assigned(self):
        response = self.job.molecules(filters=[{"type": "boolean",
                                                "value": False,
                                                "field": "assigned"}])

        self.assertEqual(response['total'], 2)

    def test_filteredon_assigned(self):
        response = self.job.molecules(filters=[{"type": "boolean",
                                                "value": True,
                                                "field": "assigned"}])

        self.assertEqual(response['total'], 0)

    def test_filteredon_reaction(self):
        filters = [{"type": "reaction",
                    "product": 3,
                    "name": u"esterase",
                    "field": "reactionsequence",
                    }]
        response = self.job.molecules(filters=filters)

        self.assertEqual(response['total'], 0)

    def test_sort_probmet(self):
        response = self.job.molecules(sorts=[{"property": "refscore",
                                              "direction": "DESC"},
                                             {"property": "molid",
                                              "direction": "ASC"}])

        self.assertEqual(response['total'], 2)

    def test_sort_nrscans(self):
        response = self.job.molecules(sorts=[{"property": "nhits",
                                              "direction": "DESC"}])

        self.assertEqual(response['total'], 2)

    def test_sort_assigned(self):
        response = self.job.molecules(sorts=[{"property": "assigned",
                                              "direction": "DESC"}])

        self.assertEqual(response['total'], 2)

    def test_sort_score(self):
        sorts = [{"property": "score", "direction": "DESC"}]

        response = self.job.molecules(scanid=641, sorts=sorts)

        self.assertEqual(response['total'], 1)

    def test_sort_score_without_scan(self):
        with self.assertRaises(ScanRequiredError):
            self.job.molecules(sorts=[{"property": "score",
                                       "direction": "DESC"}])

    def test_sort_deltappm(self):
        sorts = [{"property": "deltappm", "direction": "DESC"}]

        response = self.job.molecules(scanid=641, sorts=sorts)

        self.assertEqual(response['total'], 1)

    def test_sort_deltappm_without_scan(self):
        sorts = [{"property": "deltappm", "direction": "DESC"}]

        with self.assertRaises(ScanRequiredError):
            self.job.molecules(sorts=sorts)

    def test_sort_mz(self):
        sorts = [{"property": "mz", "direction": "DESC"}]

        response = self.job.molecules(scanid=641, sorts=sorts)

        self.assertEqual(response['total'], 1)

    def test_sort_mz_without_scan(self):
        sorts = [{"property": "mz", "direction": "DESC"}]

        with self.assertRaises(ScanRequiredError):
            self.job.molecules(sorts=sorts)


class JobDbMoleculesReactionFilterTestCase(JobDbTestCaseAbstract):

    def setUp(self):
        JobDbTestCaseAbstract.setUp(self)
        self.query = self.session.query(Molecule.molid)

    def test_reaction_reactants(self):
        afilter = {"type": "reaction",
                   "product": 3,
                   "field": "reactionsequence"}

        fq = self.job.reaction_filter(self.query, afilter)

        expected = 'SELECT molecules.molid AS molecules_molid \n'
        expected += 'FROM molecules '
        expected += 'JOIN reactions ON molecules.molid = reactions.reactant \n'
        expected += 'WHERE reactions.product = :product_1'
        self.assertEqual(str(fq), expected)

    def test_reaction_products(self):
        afilter = {"type": "reaction",
                   "reactant": 3,
                   "field": "reactionsequence"}

        fq = self.job.reaction_filter(self.query, afilter)

        expected = 'SELECT molecules.molid AS molecules_molid \n'
        expected += 'FROM molecules '
        expected += 'JOIN reactions ON molecules.molid = reactions.product \n'
        expected += 'WHERE reactions.reactant = :reactant_1'
        self.assertEqual(str(fq), expected)

    def test_reaction_productsofname(self):
        afilter = {"type": "reaction",
                   "reactant": 3,
                   "name": "esterase",
                   "field": "reactionsequence"}

        fq = self.job.reaction_filter(self.query, afilter)

        expected = 'SELECT molecules.molid AS molecules_molid \n'
        expected += 'FROM molecules '
        expected += 'JOIN reactions ON molecules.molid = reactions.product \n'
        expected += 'WHERE reactions.reactant = :reactant_1 '
        expected += 'AND reactions.name = :name_1'
        self.assertEqual(str(fq), expected)

    def test_reaction_reactantsofname(self):
        afilter = {"type": "reaction",
                   "product": 3,
                   "name": "esterase",
                   "field": "reactionsequence"}

        fq = self.job.reaction_filter(self.query, afilter)

        expected = 'SELECT molecules.molid AS molecules_molid \n'
        expected += 'FROM molecules '
        expected += 'JOIN reactions ON molecules.molid = reactions.reactant \n'
        expected += 'WHERE reactions.product = :product_1 '
        expected += 'AND reactions.name = :name_1'
        self.assertEqual(str(fq), expected)

    def test_reaction_reactantandproduct(self):
        afilter = {"type": "reaction",
                   "product": 3,
                   "reactant": 3,
                   "name": "esterase",
                   "field": "reactionsequence"}

        with self.assertRaises(TypeError):
            self.job.reaction_filter(self.query, afilter)

    def test_reaction_onlyname(self):
        afilter = {"type": "reaction",
                   "name": "esterase",
                   "field": "reactionsequence"}

        with self.assertRaises(TypeError):
            self.job.reaction_filter(self.query, afilter)

    def test_none(self):
        afilter = {"type": "reaction",
                   "field": "reactionsequence"}

        with self.assertRaises(TypeError):
            self.job.reaction_filter(self.query, afilter)


class JobDbMolecules2csvTestCase(JobDbTestCaseAbstract):

    def test_it(self):
        csvfile = self.job.molecules2csv(self.job.molecules()['rows'])
        import csv
        import StringIO
        expected_csvfile = StringIO.StringIO()
        cols = ['name', 'smiles', 'refscore', 'reactionsequence',
                'nhits', 'formula', 'mim', 'predicted', 'logp', 'reference']
        csvwriter = csv.DictWriter(expected_csvfile, cols)
        csvwriter.writeheader()
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += '?cid=152432">CID: 152432</a>'
        rs1 = '{"reactantof": {"esterase": {"nr": 2, "nrp": 1}}}'
        rs2 = '{"productof": {"theogallin": {"nr": 1, "nrp": 0}}}'
        csvwriter.writerow({'name': 'pyrocatechol',
                            'smiles': 'C1=CC=C(C(=C1)O)O',
                            'refscore': 1.0,
                            'reactionsequence': rs1,
                            'nhits': 1,
                            'formula': 'C6H6O2',
                            'predicted': False,
                            'mim': 110.03677,
                            'logp': 1.231,
                            'reference': url1,
                            })
        csvwriter.writerow({'name': 'dihydroxyphenyl-valerolactone',
                            'smiles': 'O=C1CCC(Cc2ccc(O)c(O)c2)O1',
                            'refscore': 1.0,
                            'reactionsequence': rs2,
                            'nhits': 1,
                            'formula': 'C11H12O4',
                            'predicted': False,
                            'mim': 208.07355,
                            'logp': 2.763,
                            'reference': url2,
                            })
        self.assertMultiLineEqual(csvfile.getvalue(),
                                  expected_csvfile.getvalue())

    def test_with_score(self):
        mets = self.job.molecules(scanid=641)['rows']
        csvfile = self.job.molecules2csv(mets)
        import csv
        import StringIO
        expected_csvfile = StringIO.StringIO()
        cols = ['name', 'smiles', 'refscore', 'reactionsequence',
                'nhits', 'formula', 'mim', 'predicted', 'logp', 'reference',
                'score']
        csvwriter = csv.DictWriter(expected_csvfile, cols)
        csvwriter.writeheader()
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        rs1 = '{"reactantof": {"esterase": {"nr": 2, "nrp": 1}}}'
        csvwriter.writerow({'name': 'pyrocatechol',
                            'smiles': 'C1=CC=C(C(=C1)O)O',
                            'refscore': 1.0,
                            'reactionsequence': rs1,
                            'nhits': 1,
                            'formula': 'C6H6O2',
                            'predicted': False,
                            'score': 200.0,
                            'mim': 110.03677,
                            'logp': 1.231,
                            'reference': url1
                            })
        self.assertMultiLineEqual(csvfile.getvalue(),
                                  expected_csvfile.getvalue())

    def test_some_columns(self):
        cols = ['name', 'mim']
        mets = self.job.molecules(scanid=641)['rows']

        response = self.job.molecules2csv(mets, cols=cols)

        self.assertEquals(
            response.getvalue(),
            'name,mim\r\n' +
            'pyrocatechol,110.03677\r\n'
        )


class JobMolecules2sdfTestCase(JobDbTestCaseAbstract):

    def test_it(self):
        sdffile = self.job.molecules2sdf(self.job.molecules()['rows'])
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        url2 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url2 += '?cid=152432">CID: 152432</a>'

        expected_sdf = """Molfile> <name>
pyrocatechol

> <smiles>
C1=CC=C(C(=C1)O)O

> <refscore>
1.0

> <reactionsequence>
{{"reactantof": {{"esterase": {{"nr": 2, "nrp": 1}}}}}}

> <nhits>
1

> <formula>
C6H6O2

> <mim>
110.03677

> <logp>
1.231

> <reference>
{url1}

$$$$
Molfile of dihydroxyphenyl-valerolactone> <name>
dihydroxyphenyl-valerolactone

> <smiles>
O=C1CCC(Cc2ccc(O)c(O)c2)O1

> <refscore>
1.0

> <reactionsequence>
{{"productof": {{"theogallin": {{"nr": 1, "nrp": 0}}}}}}

> <nhits>
1

> <formula>
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
        mets = self.job.molecules(scanid=641)['rows']
        sdffile = self.job.molecules2sdf(mets)
        url1 = '<a href="http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi'
        url1 += '?cid=289">CID: 289</a>'
        expected_sdf = """Molfile> <name>
pyrocatechol

> <smiles>
C1=CC=C(C(=C1)O)O

> <refscore>
1.0

> <reactionsequence>
{{"reactantof": {{"esterase": {{"nr": 2, "nrp": 1}}}}}}

> <nhits>
1

> <formula>
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
        cols = ['name', 'mim']

        mets = self.job.molecules(scanid=641)['rows']
        sdffile = self.job.molecules2sdf(mets, cols=cols)

        expected_sdf = """Molfile> <name>
pyrocatechol

> <mim>
110.03677

$$$$
"""

        self.assertMultiLineEqual(sdffile, expected_sdf)


class JobScansWithMoleculesTestCase(JobDbTestCaseAbstract):

    def test_molid(self):
        response = self.job.scansWithMolecules(molid=72)
        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_all(self):
        response = self.job.scansWithMolecules()
        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_formula(self):
        filters = [{"type": "string", "value": u"C6", "field": "formula"}]

        response = self.job.scansWithMolecules(filters=filters)

        self.assertEqual(response, [{
            'rt': 933.317,
            'id': 641
        }])

    def test_nrscans(self):
        filters = [{"type": "numeric", "value": "1",
                    "comparison": "eq", "field": "nhits"}]

        response = self.job.scansWithMolecules(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_score(self):
        filters = [{"type": "numeric", "value": "200",
                    "comparison": "eq", "field": "score"}]

        response = self.job.scansWithMolecules(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317}
        ])

    def test_filteredon_not_assigned(self):
        filters = [{"type": "boolean", "value": False, "field": "assigned"}]

        response = self.job.scansWithMolecules(filters=filters)

        self.assertEqual(response, [
            {'id': 641, 'rt': 933.317},
            {'id': 870, 'rt': 1254.15}
        ])

    def test_filteredon_assigned(self):
        filters = [{"type": "boolean", "value": True, "field": "assigned"}]

        response = self.job.scansWithMolecules(filters=filters)

        self.assertEqual(response, [])

    def test_filteredon_nrscansgt_and_formula(self):
        filters = [{"type": "numeric", "comparison": "lt",
                    "value": 2, "field": "nhits"},
                   {"type": "string", "value": u"C6", "field": "formula"}]

        response = self.job.scansWithMolecules(filters=filters)

        self.assertEqual(len(response), 1)

    def test_filteron_deltappm(self):
        filters = [{"type": "numeric", "comparison": "gt",
                    "value": 0, "field": "nhits"},
                   {"type": "numeric", "value": "200",
                    "comparison": "lt", "field": "deltappm"}]

        response = self.job.scansWithMolecules(filters=filters)

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
                     'assigned_molid': None},
                    {'intensity': 807576.625,
                     'mz': 305.033508300781,
                     'assigned_molid': None}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': {'id': None, 'mz': None},
                'fragments': [{'mz': 109.0295639038086}],
            }
        )

    def test_withmslevel(self):
        self.assertEqual(
            self.job.mspectra(641, 1),
            {
                'peaks': [
                    {'intensity': 345608.65625,
                     'mz': 109.0295639038086,
                     'assigned_molid': None},
                    {'intensity': 807576.625,
                     'mz': 305.033508300781,
                     'assigned_molid': None}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': {'id': None, 'mz': None},
                'fragments': [{'mz': 109.0295639038086}],
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
                 'assigned_molid': None},
                {'intensity': 279010.28125,
                 'mz': 163.076232910156,
                 'assigned_molid': None}
            ],
            'cutoff': 13950500.0,
            'mslevel': 2,
            'precursor': {'id': 870, 'mz': 207.0663147}
        })

    def test_withassigned_met2peak(self):
        self.job.assign_molecule2peak(641, 109.0295639038086, 72)
        self.assertEqual(
            self.job.mspectra(641),
            {
                'peaks': [
                    {'intensity': 345608.65625,
                     'mz': 109.0295639038086,
                     'assigned_molid': 72},
                    {'intensity': 807576.625,
                     'mz': 305.033508300781,
                     'assigned_molid': None}
                ],
                'cutoff': 200000.0,
                'mslevel': 1,
                'precursor': {'id': None, 'mz': None},
                'fragments': [{'mz': 109.0295639038086}],
            }
        )


class JobFragmentsTestCase(JobDbTestCaseAbstract):

    def test_moleculewithoutfragments(self):
        self.maxDiff = None
        response = self.job.fragments(molid=72, scanid=641, node='root')
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
                'molid': 72,
                'mol': u'Molfile',
                'formula': u'C5H4',
                'mslevel': 1,
                'mz': 109.0295639038086,
                'scanid': 641,
                'score': 200.0,
                'isAssigned': False,
            }], 'expanded': True
        })

    def test_moleculewithfragments(self):
        self.maxDiff = None
        response = self.job.fragments(molid=352, scanid=870, node='root')
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
                    'molid': 352,
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
                    'molid': 352,
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
                'molid': 352,
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
        response = self.job.fragments(molid=352, scanid=870, node=1709)
        self.assertEqual(response, [{
            'atoms': "4,5,6,7,8,9,11,13,14",
            'deltah': 3,
            'expanded': True,
            'fragid': 1710,
            'leaf': True,
            'mass': 116.0626002568,
            'molid': 352,
            'mol': "Molfile of dihydroxyphenyl-valerolactone",
            'mslevel': 3,
            'mz': 119.08654022216797,
            'scanid': 872,
            'score': 4,
            'deltappm': 5.0781684060061766e-08,
            'formula': 'C3H5O3',
        }])

    def test_moleculewithassignedpeak(self):
        self.job.assign_molecule2peak(641, 109.0295639038086, 72)
        response = self.job.fragments(molid=72, scanid=641, node='root')
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
                'molid': 72,
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
            self.job.fragments(molid=70002, scanid=641, node='root')


class JobWithAllPeaksTestCase(unittest.TestCase):

    def setUp(self):
        self.job = JobDb(initTestingDB(dataset='useallpeaks'))

    def test_default(self):
        response = self.job.molecules()
        self.assertEquals(
            response,
            {
                'total': 1,
                'page': 1,
                'rows': [{
                    'molid': 12,
                    'predicted': True,
                    'mol': u'Molfile',
                    'formula': u'C11H12O7S',
                    'nhits': 1,
                    'name': u'5-(3,4)-dihydroxyphenyl-g-valerolactone (F)',
                    'refscore': 0.119004,
                    'reactionsequence': [u'sulfation_(aromatic_hydroxyl)'],
                    'smiles': u'Oc1ccc(CC2OC(=O)CC2)cc1OS(O)(=O)=O',
                    'inchikey14': u'YAXFVDUJDAQPTJ',
                    'mim': 288.0303734299, 'logp':1.9027,
                    'assigned': False,
                    'reference': u''
                }]
            }
        )

    def test_lvl1fragments(self):
        response = self.job.fragments(molid=12, scanid=1, node='root')
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
                'molid': 12,
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
                'molid': 12,
                'mol': u'Molfile',
                'mslevel': 1,
                'mass': 288.0303734299,
                'scanid': 1,
                'score': 0.5,
                'isAssigned': False,
                'formula': "C4H6O2",
                'children': [{
                    'fragid': 19,
                    'molid': 12,
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
                    'molid': 12,
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
