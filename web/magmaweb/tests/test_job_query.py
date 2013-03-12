import unittest
from magmaweb.job import JobQuery

class JobQueryTestCase(unittest.TestCase):
    def setUp(self):
        self.jobdir = '/somedir'
        self.jobquery = JobQuery(directory=self.jobdir)

    def test_eq(self):
        self.assertEqual(JobQuery(directory=self.jobdir),
                         self.jobquery)

    def test_eq_dir(self):
        self.assertNotEqual(JobQuery(directory='/otherdir'),
                            self.jobquery)

    def test_eq_script(self):
        job1 = JobQuery(directory=self.jobdir, script='b')
        self.assertNotEqual(job1, self.jobquery)

    def test_eq_prestaged(self):
        job1 = JobQuery(directory=self.jobdir, prestaged=[1])
        self.assertNotEqual(job1, self.jobquery)

    def test_repr(self):
        jq = JobQuery('y', script='z', prestaged=[123],
                      status_callback_url='foo')
        s = "JobQuery('y', script='z', prestaged=[123], "
        s += "status_callback_url='foo')"
        self.assertEqual(jq.__repr__(), s)

    def test_escape_single_quote(self):
        jq = JobQuery('/y')
        self.assertEquals(jq.escape("'"), '&#39;')

    def test_defaults(self):
        expected = dict(n_reaction_steps=2,
                        metabolism_types=['phase1', 'phase2'],
                        ionisation_mode=1,
                        ms_intensity_cutoff=1000000.0,
                        msms_intensity_cutoff=10,
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
        import tempfile
        self.jobdir = tempfile.mkdtemp()
        self.jobquery = JobQuery(directory=self.jobdir, script='')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.jobdir)

    def fetch_file(self, filename):
        import os.path
        return file(os.path.join(self.jobdir, filename)).read()


class JobQueryAddStructuresTestCase(JobQueryActionTestCase):
    def test_structures_as_string(self):
        params = {'structure_format': 'smiles', 'structures': 'CCO Ethanol'}
        query = self.jobquery.add_structures(params)

        sf = 'structures.dat'
        script = "{magma} add_structures -t 'smiles' structures.dat {db}\n"
        expected_query = JobQuery(directory=self.jobdir,
                                  prestaged=[sf],
                                  script=script
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

        script = "{magma} add_structures -t 'smiles' structures.dat {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['structures.dat'],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('structures.dat'))

    def test_with_metabolize(self):
        from webob.multidict import MultiDict
        params = MultiDict(structure_format='smiles',
                           structures='CCO Ethanol',
                           metabolize='on',
                           n_reaction_steps=2,
                           metabolism_types='phase1'
                           )
        params.add('metabolism_types', 'phase2')

        query = self.jobquery.add_structures(params)

        sf = 'structures.dat'
        script = "{magma} add_structures -t 'smiles' structures.dat {db} |"
        script += "{magma} metabolize --n_reaction_steps '2' "
        script += "-m 'phase1,phase2' -j - {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [sf],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

    def test_with_annotate(self):
        params = {'structure_format': 'smiles',
                  'structures': 'CCO Ethanol',
                  'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  }
        query = self.jobquery.add_structures(params, True)

        sf = 'structures.dat'
        script = "{magma} add_structures -t 'smiles' structures.dat {db} |"
        script += "{magma} annotate -p '5.0' -q '0.001' -c '200000.0'"
        script += " -d '10.0'"
        script += " -i '1'"
        script += " -b '4' --precursor_mz_precision '0.005'"
        script += " --max_water_losses '1' -j - --fast {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [sf],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

    def test_with_metabolize_and_annotate(self):
        from webob.multidict import MultiDict
        params = MultiDict(structure_format='smiles',
                           structures='CCO Ethanol',
                           metabolize='on',
                           n_reaction_steps=2,
                           metabolism_types='phase2',
                           precursor_mz_precision=0.005,
                           mz_precision=5.0,
                           mz_precision_abs=0.001,
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=10,
                           ionisation_mode=1,
                           max_broken_bonds=4,
                           max_water_losses=1,
                           )
        query = self.jobquery.add_structures(params, True)

        sf = 'structures.dat'
        script = "{magma} add_structures -t 'smiles' structures.dat {db} |"
        script += "{magma} metabolize --n_reaction_steps '2' "
        script += "-m 'phase2' -j - {db} |"
        script += "{magma} annotate -p '5.0' -q '0.001' -c '200000.0'"
        script += " -d '10.0'"
        script += " -i '1'"
        script += " -b '4' --precursor_mz_precision '0.005'"
        script += " --max_water_losses '1' -j - --fast {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [sf],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

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
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -l '3' -a '1000.0' ms_data.dat {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat'],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_ms_data_as_string(self):
        params = {'ms_data_format': 'mzxml',
                  'ms_data': 'foo',
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -l '3' -a '1000.0' ms_data.dat {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat'],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_without_ms_data(self):
        params = {'ms_data_format': 'mzxml',
                  'max_ms_level': 3,
                  'abs_peak_cutoff': 1000,
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
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  }

        query = self.jobquery.add_ms_data(params, True)

        script = "{magma} read_ms_data --ms_data_format 'mzxml' "
        script += "-l '3' -a '1000.0' ms_data.dat {db}\n"
        script += "{magma} annotate -p '5.0' -q '0.001' -c '200000.0'"
        script += " -d '10.0'"
        script += " -i '1' -b '4' --precursor_mz_precision '0.005'"
        script += " --max_water_losses '1' --fast {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat'],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_with_mass_tree_format(self):
        params = {'ms_data_format': 'mass_tree',
                  'ms_data': 'foo',
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mass_tree'"
        script += " -l '0' -a '0.0' ms_data.dat {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat'],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

    def test_with_form_tree_pos_format(self):
        params = {'ms_data_format': 'form_tree',
                  'ionisation_mode': "1",
                  'ms_data': 'foo',
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'form_tree_pos'"
        script += " -l '0' -a '0.0' ms_data.dat {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat'],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

    def test_with_form_tree_neg_format(self):
        params = {'ms_data_format': 'form_tree',
                  'ionisation_mode': "-1",
                  'ms_data': 'foo',
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'form_tree_neg'"
        script += " -l '0' -a '0.0' ms_data.dat {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat'],
                                     'script': script
                                     })
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
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        script += " -l '3' -a '1000.0' --scan '5' ms_data.dat {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat'],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_non_mzxml_with_scan(self):
        params = {'ms_data_format': 'mass_tree',
                  'ms_data': 'foo',
                  'scan': 5,
                  }

        query = self.jobquery.add_ms_data(params)

        script = "{magma} read_ms_data --ms_data_format 'mass_tree'"
        script += " -l '0' -a '0.0' ms_data.dat {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat'],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

class JobQueryMetabolizeTestCase(JobQueryActionTestCase):

    def test_it(self):
        from webob.multidict import MultiDict
        params = MultiDict(n_reaction_steps=2, metabolism_types='phase1')

        query = self.jobquery.metabolize(params)

        script = "{magma} metabolize --n_reaction_steps '2' -m 'phase1' {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

    def test_with_annotate(self):
        from webob.multidict import MultiDict
        params = MultiDict([('n_reaction_steps', 2),
                            ('metabolism_types', 'phase1'),
                            ('metabolism_types', 'phase2'),
                            ('precursor_mz_precision', 0.005),
                            ('mz_precision', 5.0),
                            ('mz_precision_abs', 0.001),
                            ('ms_intensity_cutoff', 200000),
                            ('msms_intensity_cutoff', 10),
                            ('ionisation_mode', 1),
                            ('max_broken_bonds', 4),
                            ('max_water_losses', 1),
                            ])
        query = self.jobquery.metabolize(params, True)

        script = "{magma} metabolize --n_reaction_steps '2' "
        script += "-m 'phase1,phase2' {db} |"
        script += "{magma} annotate -p '5.0' -q '0.001'"
        script += " -c '200000.0' -d '10.0'"
        script += " -i '1' -b '4' --precursor_mz_precision '0.005' "
        script += "--max_water_losses '1' -j - "
        script += "--fast {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)


class JobQueryMetabolizeOneTestCase(JobQueryActionTestCase):

    def test_it(self):
        from webob.multidict import MultiDict
        params = MultiDict(metid=123,
                           n_reaction_steps=2,
                           metabolism_types='phase1'
                           )
        params.add('metabolism_types', 'phase2')

        query = self.jobquery.metabolize_one(params)

        script = "echo '123' | {magma} metabolize -j - --n_reaction_steps '2'"
        script += " -m 'phase1,phase2' {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

    def test_with_annotate(self):
        from webob.multidict import MultiDict
        params = MultiDict(metid=123,
                           n_reaction_steps=2,
                           metabolism_types='phase1',
                           precursor_mz_precision=0.005,
                           mz_precision=5.0,
                           mz_precision_abs=0.001,
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=10,
                           ionisation_mode=1,
                           max_broken_bonds=4,
                           max_water_losses=1,
                           )

        query = self.jobquery.metabolize_one(params, True)

        script = "echo '123' | {magma} metabolize -j - --n_reaction_steps '2' "
        script += "-m 'phase1' {db} |"
        script += "{magma} annotate -p '5.0' -q '0.001' "
        script += "-c '200000.0' -d '10.0' "
        script += "-i '1' -b '4' --precursor_mz_precision '0.005'"
        script += " --max_water_losses '1' -j - "
        script += "--fast {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)


class JobQueryAnnotateTestCase(JobQueryActionTestCase):

    def test_all_params(self):
        params = {'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  }

        query = self.jobquery.annotate(params)

        script = "{magma} annotate -p '5.0' -q '0.001' -c '200000.0' -d '10.0'"
        script += " -i '1'"
        script += " -b '4' --precursor_mz_precision '0.005'"
        script += " --max_water_losses '1' --fast {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

    def test_missing_params(self):
        params = {'ionisation_mode': 1,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  }

        query = self.jobquery.annotate(params)

        script = "{magma} annotate -p '0.0' -q '0.0' -c '0.0' -d '0.0'"
        script += " -i '1'"
        script += " -b '4' --precursor_mz_precision '0.0'"
        script += " --max_water_losses '1' --fast {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [],
                                     'script': script
                                     })
        self.assertEqual(query, expected_query)

    def test_with_structure_database_without_location(self):
        params = {'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  'structure_database': 'pubchem',
                  'min_refscore': 1,
                  'max_mz': 1200,
                  }

        from colander import Invalid
        with self.assertRaises(Invalid) as e:
            self.jobquery.annotate(params)

        msg = 'Unable to locate structure database'
        self.assertEquals(e.exception.msg, msg)

    def test_with_structure_database(self):
        params = {'precursor_mz_precision': 0.005,
                  'mz_precision': 5.0,
                  'mz_precision_abs': 0.001,
                  'ms_intensity_cutoff': 200000,
                  'msms_intensity_cutoff': 10,
                  'ionisation_mode': 1,
                  'max_broken_bonds': 4,
                  'max_water_losses': 1,
                  'structure_database': 'pubchem',
                  'min_refscore': 1,
                  'max_mz': 1200,
                  }

        structure_db_location = 'data/pubchem.db'

        query = self.jobquery.annotate(params, False, structure_db_location)

        script = "{magma} annotate -p '5.0' -q '0.001' -c '200000.0' -d '10.0'"
        script += " -i '1' -b '4' --precursor_mz_precision '0.005'"
        script += " --max_water_losses '1'"
        script += " --structure_database 'pubchem'"
        script += " --db_options 'data/pubchem.db,1200,False,1'"
        script += " --fast {db}\n"
        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': [],
                                     'script': script
                                     })
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
        from webob.multidict import MultiDict
        params = MultiDict(n_reaction_steps=2,
                           ionisation_mode=1,
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=10,
                           abs_peak_cutoff=1000,
                           precursor_mz_precision=0.005,
                           max_broken_bonds=4,
                           max_water_losses=1,
                           mz_precision=5.0,
                           mz_precision_abs=0.001,
                           metabolize='on',
                           metabolism_types='phase1',
                           max_ms_level=3,
                           structures='C1CCCC1 comp1',
                           ms_data_file=msfield,
                           structure_format='smiles',
                           ms_data_format='mzxml',
                           )
        params.add('metabolism_types', 'phase2')

        query = self.jobquery.allinone(params)

        expected_script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        expected_script += " -l '3' -a '1000.0' ms_data.dat {db}\n"

        expected_script += "{magma} add_structures -t 'smiles'"
        expected_script += " structures.dat {db}\n"

        expected_script += "{magma} metabolize --n_reaction_steps '2'"
        expected_script += " -m 'phase1,phase2' {db}\n"

        expected_script += "{magma} annotate -p '5.0' -q '0.001' -c '200000.0'"
        expected_script += " -d '10.0'"
        expected_script += " -i '1' -b '4' --precursor_mz_precision '0.005'"
        expected_script += " --max_water_losses '1' --fast {db}\n"

        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat',
                                                   'structures.dat'],
                                     'script': expected_script
                                     })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(params['structures'],
                                  self.fetch_file('structures.dat'))
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_without_metabolize(self):
        self.maxDiff = 100000
        import tempfile
        from cgi import FieldStorage
        ms_data_file = tempfile.NamedTemporaryFile()
        ms_data_file.write('foo')
        ms_data_file.flush()
        msfield = FieldStorage()
        msfield.file = ms_data_file
        from webob.multidict import MultiDict
        params = MultiDict(n_reaction_steps=2,
                           ionisation_mode=1,
                           ms_intensity_cutoff=200000,
                           msms_intensity_cutoff=10,
                           abs_peak_cutoff=1000,
                           precursor_mz_precision=0.005,
                           max_broken_bonds=4,
                           max_water_losses=1,
                           mz_precision=5.0,
                           mz_precision_abs=0.001,
                           metabolism_types='phase1',
                           max_ms_level=3,
                           structures='C1CCCC1 comp1',
                           ms_data_file=msfield,
                           structure_format='smiles',
                           ms_data_format='mzxml',
                           structure_database='',
                           min_refscore=1,
                           max_mz=9999,
                           )
        params.add('metabolism_types', 'phase2')

        query = self.jobquery.allinone(params)

        expected_script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        expected_script += " -l '3' -a '1000.0' ms_data.dat {db}\n"

        expected_script += "{magma} add_structures -t 'smiles'"
        expected_script += " structures.dat {db}\n"

        expected_script += "{magma} annotate -p '5.0' -q '0.001' -c '200000.0'"
        expected_script += " -d '10.0'"
        expected_script += " -i '1' -b '4' --precursor_mz_precision '0.005'"
        expected_script += " --max_water_losses '1' --fast {db}\n"

        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat',
                                                   'structures.dat'],
                                     'script': expected_script
                                     })
        self.assertEqual(query, expected_query)
        self.assertMultiLineEqual(params['structures'],
                                  self.fetch_file('structures.dat'))
        self.assertMultiLineEqual('foo', self.fetch_file('ms_data.dat'))

    def test_without_molecule_and_with_structure_database(self):
        params = dict(n_reaction_steps=2,
                      ionisation_mode=1,
                      ms_intensity_cutoff=200000,
                      msms_intensity_cutoff=10,
                      abs_peak_cutoff=1000,
                      precursor_mz_precision=0.005,
                      max_broken_bonds=4,
                      max_water_losses=1,
                      mz_precision=5.0,
                      mz_precision_abs=0.001,
                      metabolism_types='phase1',
                      max_ms_level=3,
                      structures='',
                      ms_data='bla',
                      structure_format='smiles',
                      ms_data_format='mzxml',
                      structure_database='pubchem',
                      min_refscore=1,
                      max_mz=1200,
                      )

        structure_db_location = 'data/pubchem.db'

        query = self.jobquery.allinone(params, structure_db_location)

        expected_script = "{magma} read_ms_data --ms_data_format 'mzxml'"
        expected_script += " -l '3' -a '1000.0' ms_data.dat {db}\n"

        expected_script += "{magma} annotate -p '5.0' -q '0.001' -c '200000.0'"
        expected_script += " -d '10.0'"
        expected_script += " -i '1' -b '4' --precursor_mz_precision '0.005'"
        expected_script += " --max_water_losses '1'"
        expected_script += " --structure_database 'pubchem'"
        expected_script += " --db_options 'data/pubchem.db,1200,False,1'"
        expected_script += " --fast {db}\n"

        expected_query = JobQuery(**{'directory': self.jobdir,
                                     'prestaged': ['ms_data.dat'],
                                     'script': expected_script,
                                     })
        self.assertEqual(query, expected_query)

    def test_without_molecule_and_structure_database(self):
        params = dict(n_reaction_steps=2,
                      ionisation_mode=1,
                      ms_intensity_cutoff=200000,
                      msms_intensity_cutoff=10,
                      abs_peak_cutoff=1000,
                      precursor_mz_precision=0.005,
                      max_broken_bonds=4,
                      max_water_losses=1,
                      mz_precision=5.0,
                      mz_precision_abs=0.001,
                      metabolism_types='phase1',
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
