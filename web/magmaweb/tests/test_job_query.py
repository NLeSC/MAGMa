"""Tests for magmaweb.job.JobQuery"""
import unittest
import tempfile
import shutil
import os.path
from webob.multidict import MultiDict
from magmaweb.job import JobQuery


class JobQueryTestCase(unittest.TestCase):

    def setUp(self):
        self.jobdir = '/somedir'
        self.jobquery = JobQuery(self.jobdir)

    def test_eq(self):
        job1 = JobQuery(self.jobdir)
        self.assertEqual(job1, self.jobquery)

    def test_eq_dir(self):
        job1 = JobQuery('/otherdir')
        self.assertNotEqual(job1, self.jobquery)

    def test_eq_script(self):
        job1 = JobQuery(self.jobdir,
                        script='b',
                        )
        self.assertNotEqual(job1, self.jobquery)

    def test_eq_prestaged(self):
        job1 = JobQuery(self.jobdir,
                        prestaged=[1],
                        )
        self.assertNotEqual(job1, self.jobquery)

    def test_eq_status_callback_url(self):
        job1 = JobQuery(self.jobdir,
                        status_callback_url='/somewhere',
                        )
        self.assertNotEqual(job1, self.jobquery)

    def test_eq_restricted(self):
        job1 = JobQuery(self.jobdir,
                        restricted=True,
                        )
        self.assertNotEqual(job1, self.jobquery)

    def test_repr(self):
        jq = JobQuery('y', script='z', prestaged=[123],
                      status_callback_url='foo', restricted=True)
        s = "JobQuery('y', script='z', prestaged=[123], "
        s += "status_callback_url='foo',"
        s += "restricted=True)"
        self.assertEqual(jq.__repr__(), s)

    def test_escape_single_quote(self):
        jq = JobQuery('/y')
        self.assertEquals(jq.escape("'"), '&#39;')

    def test_defaults(self):
        expected = dict(scenario=[
            {'type': 'phase1', 'steps': '2'},
            {'type': 'phase2', 'steps': '1'}
        ],
            ionisation_mode=1,
            ms_data_format='mzxml',
            ms_data_area='',
            ms_intensity_cutoff=0.0,
            msms_intensity_cutoff=5,
            mz_precision=5.0,
            mz_precision_abs=0.001,
            abs_peak_cutoff=5000,
            max_ms_level=10,
            precursor_mz_precision=0.005,
            max_broken_bonds=3,
            max_water_losses=1,
        )
        self.assertDictEqual(expected, JobQuery.defaults())

    def test_defaults_example(self):
        example_tree = [
            '353.087494: 69989984 (',
            '    191.055756: 54674544 (',
            '        85.029587: 2596121,',
            '        93.034615: 1720164,',
            '        109.029442: 917026,',
            '        111.045067: 1104891 (',
            '            81.034691: 28070,',
            '            83.014069: 7618,',
            '            83.050339: 25471,',
            '            93.034599: 36300,',
            '            96.021790: 8453',
            '            ),',
            '        127.039917: 2890439 (',
            '            57.034718: 16911,',
            '            81.034706: 41459,',
            '            83.050301: 35131,',
            '            85.029533: 236887,',
            '            99.045074: 73742,',
            '            109.029404: 78094',
            '            ),',
            '        171.029587: 905226,',
            '        173.045212: 2285841 (',
            '            71.013992: 27805,',
            '            93.034569: 393710,',
            '            111.008629: 26219,',
            '            111.045029: 339595,',
            '            137.024292: 27668,',
            '            155.034653: 145773',
            '            ),',
            '        191.055725: 17000514',
            '        ),',
            '    353.087097: 4146696',
            '    )'
        ]
        expected = dict(ms_data="\n".join(example_tree),
                        ms_data_format='mass_tree',
                        ionisation_mode=-1,
                        )
        self.assertDictEqual(expected, JobQuery.defaults('example'))

    def test_defaults_example2(self):
        example_tree = [
            'C16H17O9: 69989984 (',
            '    C7H11O6: 54674544 (',
            '        C4H5O2: 2596121,',
            '        C6H5O: 1720164,',
            '        C6H5O2: 917026,',
            '        C6H7O2: 1104891 (',
            '            C5H5O: 28070,',
            '            C4H3O2: 7618,',
            '            C5H7O: 25471,',
            '            C6H5O: 36300,',
            '            C5H4O2: 8453',
            '            ),',
            '        C6H7O3: 2890439 (',
            '            C3H5O: 16911,',
            '            C5H5O: 41459,',
            '            C5H7O: 35131,',
            '            C4H5O2: 236887,',
            '            C5H7O2: 73742,',
            '            C6H5O2: 78094',
            '            ),',
            '        C7H7O5: 905226,',
            '        C7H9O5: 2285841 (',
            '            C3H3O2: 27805,',
            '            C6H5O: 393710,',
            '            C5H3O3: 26219,',
            '            C6H7O2: 339595,',
            '            C7H5O3: 27668,',
            '            C7H7O4: 145773',
            '            ),',
            '        C7H11O6: 17000514',
            '        ),',
            '    C16H17O9: 4146696',
            '    )',
        ]
        expected = dict(ms_data="\n".join(example_tree),
                        ms_data_format='form_tree',
                        ionisation_mode=-1,
                        )
        self.assertDictEqual(expected, JobQuery.defaults('example2'))


class JobQueryFileTestCase(unittest.TestCase):

    def test_valid(self):
        from cgi import FieldStorage
        f = FieldStorage()
        df = JobQuery.File().deserialize(None, f)

        self.assertEquals(f, df)

    def test_null(self):
        from colander import null
        self.assertEquals(JobQuery.File().deserialize(None, null), null)

    def test_emptystring(self):
        from colander import null
        self.assertEquals(JobQuery.File().deserialize(None, ''), null)

    def test_invalid(self):
        from colander import Invalid, SchemaNode
        n = SchemaNode(JobQuery.File(), name='filefield')
        with self.assertRaises(Invalid) as e:
            n.deserialize(12345)

        self.assertDictEqual(e.exception.asdict(),
                             {'filefield': '12345 is not a cgi.FieldStorage'})

    def test_serialize(self):
        self.assertEquals(JobQuery.File().serialize(None, 12345), 12345)


class JobQueryActionTestCase(unittest.TestCase):

    def setUp(self):
        self.jobdir = tempfile.mkdtemp()
        self.jobquery = JobQuery(directory=self.jobdir,
                                 script='',
                                 status_callback_url='/',
                                 )

    def tearDown(self):
        shutil.rmtree(self.jobdir)

    def fetch_file(self, filename):
        return file(os.path.join(self.jobdir, filename)).read()


class JobQueryAddStructuresTestCase(JobQueryActionTestCase):

    def test_structures_as_string(self):
        params = {'structure_format': 'smiles', 'structures': 'CCO Ethanol'}
        query = self.jobquery.add_structures(params)

        sf = 'structures.dat'
        script = "{magma} add_structures -p -t 'smiles' structures.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[sf],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(params['structures'], self.fetch_file(sf))

    def test_structures_as_file(self):
        import tempfile
        from cgi import FieldStorage
        sfile = tempfile.TemporaryFile()
        sfile.write('foo')
        sfile.flush()
        sfield = FieldStorage()
        sfield.file = sfile
        params = {'structure_format': 'smiles', 'structures_file': sfield}

        query = self.jobquery.add_structures(params)

        script = "{magma} add_structures -p -t 'smiles' structures.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['structures.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('structures.dat'))

    def test_with_metabolize(self):

        params = MultiDict(structure_format='smiles',
                           structures='CCO Ethanol',
                           metabolize='on',
                           scenario=[
                               {'type': 'phase1', 'steps': '2'},
                               {'type': 'phase2', 'steps': '1'}
                           ],
                           )

        query = self.jobquery.add_structures(params)

        sf = 'structures.dat'
        scen = 'scenario.csv'
        script = "{magma} add_structures -p -t 'smiles' structures.dat {db} |"
        script += "{magma} metabolize -p --scenario scenario.csv"
        script += " --call_back_url '/' -j - {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[sf, scen],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(
            'CCO Ethanol', self.fetch_file('structures.dat'))
        self.assertMultiLineEqual(
            'phase1,2\nphase2,1\n', self.fetch_file('scenario.csv'))

    def test_with_annotate(self):
        params = {'structure_format': 'smiles',
                  'structures': 'CCO Ethanol',
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  }
        query = self.jobquery.add_structures(params, True)

        sf = 'structures.dat'
        script = "{magma} add_structures -p -t 'smiles' structures.dat {db} |"
        script += "{magma} annotate -c '200000.0'"
        script += " -d '10.0' -b '4'"
        script += " --max_water_losses '1' --call_back_url '/'"
        script += " -j - --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[sf],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(
            'CCO Ethanol', self.fetch_file('structures.dat'))

    def test_with_metabolize_and_annotate(self):

        params = MultiDict(structure_format='smiles',
                           structures='CCO Ethanol',
                           metabolize='on',
                           scenario=[{'type': 'phase1', 'steps': '2'},
                                     {'type': 'phase2', 'steps': '1'}],
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=10,
                           max_broken_bonds=4,
                           max_water_losses=1,
                           )
        query = self.jobquery.add_structures(params, True)

        sf = 'structures.dat'
        script = "{magma} add_structures -p -t 'smiles'"
        script += " structures.dat {db} |"
        script += "{magma} metabolize -p --scenario scenario.csv"
        script += " --call_back_url '/' -j - {db} |"
        script += "{magma} annotate -c '200000.0'"
        script += " -d '10.0' -b '4'"
        script += " --max_water_losses '1' --call_back_url '/'"
        script += " -j - --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[sf, 'scenario.csv'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(
            'CCO Ethanol', self.fetch_file('structures.dat'))
        self.assertMultiLineEqual(
            'phase1,2\nphase2,1\n', self.fetch_file('scenario.csv'))

    def test_without_structures(self):
        params = {'structure_format': 'smiles',
                  'structures_file': '',
                  'structures': ''
                  }
        from colander import Invalid
        with self.assertRaises(Invalid) as e:
            self.jobquery.add_structures(params, False)

        s = 'Either structures or structures_file must be set'
        sf = 'Either structures or structures_file must be set'
        expected = {'structures': s, 'structures_file': sf}
        self.assertDictEqual(e.exception.asdict(), expected)

    def test_with_structure_as_string_and_file(self):
        import tempfile
        from cgi import FieldStorage
        sfile = tempfile.TemporaryFile()
        sfile.write('foo')
        sfile.flush()
        sfield = FieldStorage()
        sfield.file = sfile
        params = {'structure_format': 'smiles',
                  'structures_file': sfield,
                  'structures': 'bar'
                  }

        self.jobquery.add_structures(params)
        # File is kept and string is ignored
        self.assertMultiLineEqual('foo', self.fetch_file('structures.dat'))


class JobQueryAddMSDataTestCase(JobQueryActionTestCase):

    def test_ms_data_as_file(self):
        import tempfile
        from cgi import FieldStorage
        msfile = tempfile.TemporaryFile()
        msfile.write('foo')
        msfile.flush()
        msfield = FieldStorage()
        msfield.file = msfile
        params = {'ms_data_format': 'mzxml',
                  'ms_data_file': msfield,
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -i '1' -m '3' -a '1000.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_ms_data_as_string(self):
        params = {'ms_data_format': 'mzxml',
                  'ms_data': 'foo',
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -i '1' -m '3' -a '1000.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_without_ms_data(self):
        params = {'ms_data_format': 'mzxml',
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  }
        from colander import Invalid
        with self.assertRaises(Invalid) as e:
            self.jobquery.add_ms_data(params, False)

        s = 'Either ms_data or ms_data_file must be set'
        expected = {'ms_data': s, 'ms_data_file': s}
        self.assertDictEqual(e.exception.asdict(), expected)

    def test_with_ms_data_as_string_and_file(self):
        import tempfile
        from cgi import FieldStorage
        msfile = tempfile.TemporaryFile()
        msfile.write('foo')
        msfile.flush()
        msfield = FieldStorage()
        msfield.file = msfile
        params = {'ms_data_format': 'mzxml',
                  'ms_data_file': msfield,
                  'ms_data': 'bar',
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  }

        self.jobquery.add_ms_data(params)

        # File is kept and string is ignored
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_with_annotate(self):
        import tempfile
        from cgi import FieldStorage
        msfile = tempfile.TemporaryFile()
        msfile.write('foo')
        msfile.flush()
        msfield = FieldStorage()
        msfield.file = msfile
        params = {'ms_data_format': 'mzxml',
                  'ms_data_file': msfield,
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  }

        query = self.jobquery.add_ms_data(params, True)

        script = "{magma} read_ms_data --ms_data_format 'mzxml' "
        script += "-i '1' -m '3' -a '1000.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --call_back_url '/' ms_data.dat {db}\n"
        script += "{magma} annotate -c '200000.0'"
        script += " -d '10.0'"
        script += " -b '4'"
        script += " --max_water_losses '1' --call_back_url '/' --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_with_mass_tree_format(self):
        params = {'ms_data_format': 'mass_tree',
                  'ms_data': 'foo',
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mass_tree'"
        script += " -i '1' -m '0' -a '0.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)

    def test_with_form_tree_pos_format(self):
        params = {'ms_data_format': 'form_tree',
                  'ionisation_mode': "1",
                  'ms_data': 'foo',
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'form_tree_pos'"
        script += " -i '1' -m '0' -a '0.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)

    def test_with_form_tree_neg_format(self):
        params = {'ms_data_format': 'form_tree',
                  'ionisation_mode': "-1",
                  'ms_data': 'foo',
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'form_tree_neg'"
        script += " -i '-1' -m '0' -a '0.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)

    def test_with_form_tree_noion_format(self):
        params = {'ms_data_format': 'form_tree',
                  'ms_data': 'foo',
                  }
        from colander import Invalid
        with self.assertRaises(Invalid) as e:
            self.jobquery.add_ms_data(params)

        msg = 'Require ionisation_mode when ms_data_format=form_tree'
        self.assertEquals(e.exception.msg, msg)

    def test_mzxml_with_scan(self):
        params = {'ms_data_format': 'mzxml',
                  'ms_data': 'foo',
                  'scan': 5,
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -i '1' -m '3' -a '1000.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --scan '5' --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_non_mzxml_with_scan(self):
        params = {'ms_data_format': 'mass_tree',
                  'ms_data': 'foo',
                  'scan': 5,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mass_tree'"
        script += " -i '1' -m '0' -a '0.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)

    def test_restricted_without_scan_and_molupload(self):
        self.jobquery.restricted = True
        params = {'ms_data_format': 'mzxml',
                  'ms_data': 'foo',
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -i '1' -m '3' -a '1000.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --time_limit 1 --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  restricted=True,
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_restricted_without_scan_and_structdb(self):
        self.jobquery.restricted = True
        params = {'ms_data_format': 'mzxml',
                  'ms_data': 'foo',
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'structure_database': 'pubchem',
                  'ionisation_mode': 1,
                  }

        from colander import Invalid
        with self.assertRaises(Invalid) as e:
            self.jobquery.add_ms_data(params)

        msg = 'Require MS1 scan number'
        self.assertEquals(e.exception.msg, msg)

    def test_resticted(self):
        self.jobquery.restricted = True
        params = {'ms_data_format': 'mzxml',
                  'ms_data': 'foo',
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'scan': 5,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ionisation_mode': 1,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -i '1' -m '3' -a '1000.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --scan '5'"
        script += " --time_limit 1 --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  restricted=True,
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_noionmode_pos(self):
        params = {'ms_data_format': 'mzxml',
                  'ms_data': 'foo',
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -i '1' -m '3' -a '1000.0'"
        script += " -p '5.0' -q '0.001' --precursor_mz_precision '0.005'"
        script += " --call_back_url '/' ms_data.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))


class JobQueryMetabolizeTestCase(JobQueryActionTestCase):

    def test_it(self):

        params = MultiDict(scenario=[{'type': 'phase1', 'steps': '2'},
                                     {'type': 'phase2', 'steps': '1'}]
                           )

        query = self.jobquery.metabolize(params)

        script = "{magma} metabolize -p --scenario scenario.csv"
        script += " --call_back_url '/' {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['scenario.csv'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(
            'phase1,2\nphase2,1\n', self.fetch_file('scenario.csv'))

    def test_with_jsonified_scenario(self):
        scenario = '[{"steps": "2", "type": "phase1"}, '
        scenario += '{"steps": "1", "type": "phase2"}]'
        params = dict(scenario=scenario)

        query = self.jobquery.metabolize(params)

        script = "{magma} metabolize -p --scenario scenario.csv"
        script += " --call_back_url '/' {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['scenario.csv'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(
            'phase1,2\nphase2,1\n', self.fetch_file('scenario.csv'))

    def test_with_misformatjsonified_scenario(self):
        scenario = 'xxxxxx[{"steps": "2", "type": "phase1"}, '
        scenario += '{"steps": "1", "type": "phase2"}]'
        params = dict(scenario=scenario)

        with self.assertRaises(ValueError):
            self.jobquery.metabolize(params)

    def test_with_annotate(self):

        params = MultiDict([('scenario', [{'type': 'phase1', 'steps': '2'},
                                          {'type': 'phase2', 'steps': '1'}]),
                            ('ms_intensity_cutoff', 200000),
                            ('msms_intensity_cutoff', 10),
                            ('max_broken_bonds', 4),
                            ('max_water_losses', 1),
                            ])
        query = self.jobquery.metabolize(params, True)

        script = "{magma} metabolize -p --scenario scenario.csv"
        script += " --call_back_url '/' {db} |"
        script += "{magma} annotate"
        script += " -c '200000.0' -d '10.0'"
        script += " -b '4'"
        script += " --max_water_losses '1' --call_back_url '/' -j -"
        script += " --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['scenario.csv'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(
            'phase1,2\nphase2,1\n', self.fetch_file('scenario.csv'))

    def test_it_restricted(self):
        self.jobquery.restricted = True
        params = MultiDict(scenario=[{'type': 'phase1', 'steps': '2'},
                                     {'type': 'phase2', 'steps': '1'}]
                           )

        query = self.jobquery.metabolize(params)

        script = "{magma} metabolize -p --scenario scenario.csv"
        script += " --call_back_url '/' --time_limit 3 {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['scenario.csv'],
                                  script=script,
                                  status_callback_url='/',
                                  restricted=True,
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('phase1,2\nphase2,1\n',
                                  self.fetch_file('scenario.csv')
                                  )

    def test_with_structuredb(self):
        params = MultiDict(scenario=[{'type': 'phase1', 'steps': '2'},
                                     {'type': 'phase2', 'steps': '1'}],
                           structure_database='pubchem',
                           )

        from colander import Invalid
        with self.assertRaises(Invalid) as e:
            self.jobquery.metabolize(params)

        msg = 'Not allowed to metabolize structure database'
        self.assertEquals(e.exception.msg, msg)


class JobQueryMetabolizeOneTestCase(JobQueryActionTestCase):

    def test_it(self):

        params = MultiDict(molid=123,
                           scenario=[{'type': 'phase1', 'steps': '2'},
                                     {'type': 'phase2', 'steps': '1'}],
                           )

        query = self.jobquery.metabolize_one(params)

        script = "echo '123' | {magma} metabolize -p -j - "
        script += "--scenario scenario.csv --call_back_url '/' {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['scenario.csv'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(
            'phase1,2\nphase2,1\n', self.fetch_file('scenario.csv'))

    def test_with_annotate(self):

        params = MultiDict(molid=123,
                           scenario=[{'type': 'phase1', 'steps': '2'},
                                     {'type': 'phase2', 'steps': '1'}],
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=10,
                           max_broken_bonds=4,
                           max_water_losses=1,
                           )

        query = self.jobquery.metabolize_one(params, True)

        script = "echo '123' | {magma} metabolize -p -j - --scenario scenario.csv"
        script += " --call_back_url '/' {db} |"
        script += "{magma} annotate"
        script += " -c '200000.0' -d '10.0'"
        script += " -b '4'"
        script += " --max_water_losses '1' --call_back_url '/' -j -"
        script += " --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['scenario.csv'],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(
            'phase1,2\nphase2,1\n', self.fetch_file('scenario.csv'))


class JobQueryAnnotateTestCase(JobQueryActionTestCase):

    def test_all_params(self):
        params = {'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  }

        query = self.jobquery.annotate(params)

        script = "{magma} annotate -c '200000.0' -d '10.0'"
        script += " -b '4'"
        script += " --max_water_losses '1' --call_back_url '/' --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)

    def test_missing_params(self):
        params = {'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  }

        query = self.jobquery.annotate(params)

        script = "{magma} annotate -c '0.0' -d '0.0'"
        script += " -b '4'"
        script += " --max_water_losses '1' --call_back_url '/' --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)

    def test_with_structure_database(self):
        params = {'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  'structure_database': 'pubchem',
                  'min_refscore': 1,
                  'max_mz': 1200,
                  }

        query = self.jobquery.annotate(params, False)

        script = "{magma} annotate -c '200000.0' -d '10.0'"
        script += " -b '4'"
        script += " --max_water_losses '1' --call_back_url '/'"
        script += " --structure_database 'pubchem'"
        script += " --db_options ',1200,False,True,1'"
        script += " --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)

    def test_with_structure_database_restricted(self):
        self.jobquery.restricted = True
        params = {'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  'structure_database': 'pubchem',
                  'min_refscore': 1,
                  'max_mz': 1200,
                  }

        query = self.jobquery.annotate(params, False)

        script = "{magma} annotate -c '200000.0' -d '10.0'"
        script += " -b '4'"
        script += " --max_water_losses '1' --call_back_url '/'"
        script += " --structure_database 'pubchem'"
        script += " --db_options ',1200,True,True,1'"
        script += " --time_limit 3 --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[],
                                  script=script,
                                  status_callback_url='/',
                                  restricted=True,
                                  )
        self.assertEqual(query, expected_query)

    def test_with_structure_database_exclude_halo(self):
        params = {'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  'structure_database': 'pubchem',
                  'min_refscore': 1,
                  'max_mz': 1200,
                  'excl_halo': True,
                  }

        query = self.jobquery.annotate(params, False)

        script = "{magma} annotate -c '200000.0' -d '10.0'"
        script += " -b '4'"
        script += " --max_water_losses '1' --call_back_url '/'"
        script += " --structure_database 'pubchem'"
        script += " --db_options ',1200,False,False,1'"
        script += " --fast {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[],
                                  script=script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)


class JobQueryAllInOneTestCase(JobQueryActionTestCase):

    def test_with_metabolize(self):
        self.maxDiff = 100000
        import tempfile
        from cgi import FieldStorage
        ms_data_file = tempfile.NamedTemporaryFile()
        ms_data_file.write('foo')
        ms_data_file.flush()
        msfield = FieldStorage()
        msfield.file = ms_data_file

        params = MultiDict(ionisation_mode=1,
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=10,
                           abs_peak_cutoff=1000,
                           precursor_mz_precision=0.005,
                           max_broken_bonds=4,
                           max_water_losses=1,
                           mz_precision=5.0,
                           mz_precision_abs=0.001,
                           metabolize='on',
                           scenario=[{'type': 'phase1', 'steps': '2'},
                                     {'type': 'phase2', 'steps': '1'}],
                           max_ms_level=3,
                           structures='C1CCCC1 comp1',
                           ms_data_file=msfield,
                           structure_format='smiles',
                           ms_data_format='mzxml',
                           )

        query = self.jobquery.allinone(params)

        expected_script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        expected_script += " -i '1' -m '3' -a '1000.0'"
        expected_script += " -p '5.0' -q '0.001'"
        expected_script += " --precursor_mz_precision '0.005'"
        expected_script += " --call_back_url '/'"
        expected_script += " ms_data.dat {db}\n"

        expected_script += "{magma} add_structures -p -t 'smiles'"
        expected_script += " structures.dat {db}\n"

        expected_script += "{magma} metabolize -p --scenario scenario.csv"
        expected_script += " --call_back_url '/' {db}\n"

        expected_script += "{magma} annotate -c '200000.0'"
        expected_script += " -d '10.0' -b '4'"
        expected_script += " --max_water_losses '1' --call_back_url '/'"
        expected_script += " --fast {db}\n"

        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat',
                                             'structures.dat',
                                             'scenario.csv',
                                             ],
                                  script=expected_script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(params['structures'],
                                  self.fetch_file('structures.dat'))
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))
        self.assertMultiLineEqual(
            'phase1,2\nphase2,1\n', self.fetch_file('scenario.csv'))

    def test_without_metabolize(self):
        self.maxDiff = 100000
        import tempfile
        from cgi import FieldStorage
        ms_data_file = tempfile.NamedTemporaryFile()
        ms_data_file.write('foo')
        ms_data_file.flush()
        msfield = FieldStorage()
        msfield.file = ms_data_file

        params = MultiDict(ionisation_mode=1,
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=10,
                           abs_peak_cutoff=1000,
                           precursor_mz_precision=0.005,
                           max_broken_bonds=4,
                           max_water_losses=1,
                           mz_precision=5.0,
                           mz_precision_abs=0.001,
                           scenario=[{'type': 'phase1', 'steps': '2'},
                                     {'type': 'phase2', 'steps': '1'}],
                           max_ms_level=3,
                           structures='C1CCCC1 comp1',
                           ms_data_file=msfield,
                           structure_format='smiles',
                           ms_data_format='mzxml',
                           structure_database='',
                           min_refscore=1,
                           max_mz=9999,
                           )

        query = self.jobquery.allinone(params)

        expected_script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        expected_script += " -i '1' -m '3' -a '1000.0'"
        expected_script += " -p '5.0' -q '0.001'"
        expected_script += " --precursor_mz_precision '0.005'"
        expected_script += " --call_back_url '/'"
        expected_script += " ms_data.dat {db}\n"

        expected_script += "{magma} add_structures -p -t 'smiles'"
        expected_script += " structures.dat {db}\n"

        expected_script += "{magma} annotate -c '200000.0'"
        expected_script += " -d '10.0'"
        expected_script += " -b '4'"
        expected_script += " --max_water_losses '1' --call_back_url '/'"
        expected_script += " --fast {db}\n"

        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat', 'structures.dat'],
                                  script=expected_script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(params['structures'],
                                  self.fetch_file('structures.dat'))
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_without_molecule_and_with_structure_database(self):
        params = dict(ionisation_mode=1,
                      ms_intensity_cutoff=200000,
                      msms_intensity_cutoff=10,
                      abs_peak_cutoff=1000,
                      precursor_mz_precision=0.005,
                      max_broken_bonds=4,
                      max_water_losses=1,
                      mz_precision=5.0,
                      mz_precision_abs=0.001,
                      scenario=[{'type': 'phase1', 'steps': '2'},
                                {'type': 'phase2', 'steps': '1'}],
                      max_ms_level=3,
                      structures='',
                      ms_data='bla',
                      structure_format='smiles',
                      ms_data_format='mzxml',
                      structure_database='pubchem',
                      min_refscore=1,
                      max_mz=1200,
                      )

        query = self.jobquery.allinone(params)

        expected_script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        expected_script += " -i '1' -m '3' -a '1000.0'"
        expected_script += " -p '5.0' -q '0.001'"
        expected_script += " --precursor_mz_precision '0.005'"
        expected_script += " --call_back_url '/'"
        expected_script += " ms_data.dat {db}\n"

        expected_script += "{magma} annotate -c '200000.0'"
        expected_script += " -d '10.0'"
        expected_script += " -b '4'"
        expected_script += " --max_water_losses '1' --call_back_url '/'"
        expected_script += " --structure_database 'pubchem'"
        expected_script += " --db_options ',1200,False,True,1'"
        expected_script += " --fast {db}\n"

        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=['ms_data.dat'],
                                  script=expected_script,
                                  status_callback_url='/',
                                  )
        self.assertEqual(query, expected_query)

    def test_without_molecule_and_structure_database(self):
        params = dict(ionisation_mode=1,
                      ms_intensity_cutoff=200000,
                      msms_intensity_cutoff=10,
                      abs_peak_cutoff=1000,
                      precursor_mz_precision=0.005,
                      max_broken_bonds=4,
                      max_water_losses=1,
                      mz_precision=5.0,
                      mz_precision_abs=0.001,
                      scenario=[{'type': 'phase1', 'steps': '2'},
                                {'type': 'phase2', 'steps': '1'}],
                      max_ms_level=3,
                      structures='',
                      ms_data='bla',
                      structure_format='smiles',
                      ms_data_format='mzxml',
                      structure_database='',
                      min_refscore=1,
                      max_mz=9999,
                      )

        from colander import Invalid
        with self.assertRaises(Invalid) as e:
            self.jobquery.allinone(params)

        s = 'Either structures or structures_file must be set'
        sf = 'Either structures or structures_file must be set'
        sd = 'Either structures or structures_file'
        sd += ' or structure_database must be set'
        expected = {'structures': s,
                    'structures_file': sf,
                    'structure_database': sd}
        self.assertDictEqual(e.exception.asdict(), expected)
