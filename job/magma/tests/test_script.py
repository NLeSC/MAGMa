import unittest
import mock
import argparse
import magma.script

class TestMagmaCommand(unittest.TestCase):

    def setUp(self):
        self.mc = magma.script.MagmaCommand()

    def test_version(self):

        self.assertEqual(self.mc.version(), '1.0')

    def test_all_in_one(self):
        # mock MagmaSession
        ms = mock.Mock(magma.MagmaSession)
        self.mc.get_magma_session = mock.Mock(return_value=ms)
        se = mock.Mock(magma.StructureEngine)
        ms.get_structure_engine.return_value = se
        me = mock.Mock(magma.MsDataEngine)
        ms.get_ms_data_engine.return_value = me
        ae = mock.Mock(magma.AnnotateEngine)
        ms.get_annotate_engine.return_value = ae

        args = argparse.Namespace()
        args.db = ':memory:'
        args.description = 'my desc'
        args.metabolism_types = 'phase1'
        args.n_reaction_steps = 2
        args.structure_format = 'smiles'
        args.structures = [ 'CCO|ethanol', '' ]
        args.abs_peak_cutoff = 100000
        args.rel_peak_cutoff = 0.05
        args.max_ms_level = 5
        args.ms_data_format = 'mzxml'
        args.ms_data = argparse.Namespace(name='bogus.mzxml')
        args.ionisation_mode = 1
        args.skip_fragmentation = False
        args.max_broken_bonds = 4
        args.ms_intensity_cutoff = 20000
        args.msms_intensity_cutoff = 0.001
        args.mz_precision = 0.01
        args.precursor_mz_precision = 0.123
        args.use_all_peaks = True

        self.mc.all_in_one(args)

        self.mc.get_magma_session.assert_called_with(args.db, args.description)
        ms.get_structure_engine.assert_called_with(args.metabolism_types, args.n_reaction_steps)
        self.assertTrue(se.add_structure.called)
        se.metabolize_all.assert_called_with(args.metabolism_types, args.n_reaction_steps)
        ms.get_ms_data_engine.assert_called_with(abs_peak_cutoff=args.abs_peak_cutoff,
            rel_peak_cutoff=args.rel_peak_cutoff, max_ms_level=args.max_ms_level)
        me.store_mzxml_file.assert_called_with(args.ms_data.name)
        ms.get_annotate_engine.assert_called_with(
            ionisation_mode=args.ionisation_mode,
            skip_fragmentation=args.skip_fragmentation,
            max_broken_bonds=args.max_broken_bonds,
            ms_intensity_cutoff=args.ms_intensity_cutoff,
            msms_intensity_cutoff=args.msms_intensity_cutoff,
            mz_precision=args.mz_precision,
            precursor_mz_precision=args.precursor_mz_precision,
            use_all_peaks=args.use_all_peaks
        )
        ae.build_spectra.assert_called_with()
        ae.search_all_structures.assert_called_with()


