import unittest
import argparse
import magma.script
import tempfile, os
from magma.models import Base, Molecule, Scan, Peak, Fragment, Run

class TestMagmaCommand(unittest.TestCase):

    def setUp(self):
        self.mc = magma.script.MagmaCommand()

    def test_version(self):

        self.assertEqual(self.mc.version(), '1.0')

    def test_chlorogenic_acid_example_without_fast_option(self):
        treefile = tempfile.NamedTemporaryFile(delete=False)
        dbfile = tempfile.NamedTemporaryFile(delete=False)

        args = argparse.Namespace()
        args.db = dbfile.name
        args.ms_data = treefile
        args.description = 'Example'
        args.ionisation_mode = -1
        args.abs_peak_cutoff = 0
        args.mz_precision = 5
        args.mz_precision_abs = 0.001
        args.precursor_mz_precision = 0.005
        args.max_ms_level = 5
        args.ms_data_format = 'mass_tree'
        #args.ms_data = argparse.Namespace(name='bogus.mzxml')
        args.log = 'debug'
        args.call_back_url = None

        treefile.write("""353.087494: 69989984 (
    191.055756: 54674544 (
        85.029587: 2596121,
        93.034615: 1720164,
        109.029442: 917026,
        111.045067: 1104891 (
            81.034691: 28070,
            83.014069: 7618,
            83.050339: 25471,
            93.034599: 36300,
            96.021790: 8453
            ),
        127.039917: 2890439 (
            57.034718: 16911,
            81.034706: 41459,
            83.050301: 35131,
            85.029533: 236887,
            99.045074: 73742,
            109.029404: 78094
            ),
        171.029587: 905226,
        173.045212: 2285841 (
            71.013992: 27805,
            93.034569: 393710,
            111.008629: 26219,
            111.045029: 339595,
            137.024292: 27668,
            155.034653: 145773
            ),
        191.055725: 17000514
        ),
    353.087097: 4146696
    )
""")
        treefile.close()

        self.mc.read_ms_data(args)
        os.remove(treefile.name)

        args = argparse.Namespace()        
        args.db = dbfile.name
        args.description = None
        args.skip_fragmentation = False
        args.max_broken_bonds = 3
        args.max_water_losses = 1
        args.ms_intensity_cutoff = 0
        args.msms_intensity_cutoff = 0
        args.use_all_peaks = True
        args.adducts = None
        args.max_charge = 1
        args.log = 'debug'
        args.call_back_url = None
        args.scans= 'all'
        args.structure_database = 'hmdb'
        args.db_options=''
        args.molids = None
        args.ncpus = 1
        args.fast = False
        args.time_limit = None

        self.mc.annotate(args)

        ms = magma.MagmaSession(dbfile.name)
        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset(
                          {
                          'ionisation_mode': -1,
                          'abs_peak_cutoff': 0,
                          'max_ms_level': 5,
                          'mz_precision': 5.0,
                          'precursor_mz_precision': 0.005,
                          'max_broken_bonds': 3,
                          'ms_intensity_cutoff': 0
                           }, rundata.__dict__)
        scandata = ms.db_session.query(Scan).count()
        self.assertEqual(scandata,6)
        peakdata = ms.db_session.query(Peak).count()
        self.assertEqual(peakdata,28)
        moleculedata = ms.db_session.query(Molecule).count()
        self.assertGreater(moleculedata,4)
        fragmentdata = ms.db_session.query(Fragment).count()
        self.assertGreater(fragmentdata,86)

        os.remove(dbfile.name)

    def test_JWH015_example_with_fast_option(self):
        dbfile = tempfile.NamedTemporaryFile(delete=False)

        treefile = tempfile.NamedTemporaryFile(delete=False)
        args = argparse.Namespace()
        args.db = dbfile.name
        args.ms_data = treefile
        args.description = 'Example'
        args.ionisation_mode = 1
        args.abs_peak_cutoff = 0
        args.mz_precision = 10
        args.mz_precision_abs = 0.002
        args.precursor_mz_precision = 0.005
        args.max_ms_level = 5
        args.ms_data_format = 'mass_tree'
        #args.ms_data = argparse.Namespace(name='bogus.mzxml')
        args.log = 'debug'
        args.call_back_url = None

        treefile.write("""328.1696:100 (
    155.0492:1,
    127.0541:1,
    200.1069:1
    ),
286.1226:100 (
    127.0540:1800,
    130.0654:200,
    155.0489:1600,
    158.0599:800
    ),
344.1645:100 (
    127.0537:250,
    155.0494:350,
    216.1016:70,
    302.1467:40
    ),
520.1966:100 (
    155.0493:1,
    127.0537:1,
    216.1019:1,
    344.1648:1
    )
""")
        treefile.close()

        self.mc.read_ms_data(args)
        os.remove(treefile.name)

        args = argparse.Namespace()        
        args.db = dbfile.name
        args.description = None
        args.log = 'debug'
        args.pubchem_names = False
        args.structure_format = 'smiles'
        args.mass_filter = 9999
        args.structures = 'c1ccc2ccccc2c1C(=O)c3c4ccccc4n(CCC)c3C'
        self.mc.add_structures(args)

        scenariofile = tempfile.NamedTemporaryFile(delete=False)
        args = argparse.Namespace()
        args.db = dbfile.name
        args.description = None
        args.log = 'debug'
        args.pubchem_names = False
        args.call_back_url = None
        args.scenario = scenariofile.name
        args.time_limit = None
        scenariofile.write('phase1,1\nphase2,1')
        scenariofile.close()
        self.mc.metabolize(args)
        os.remove(scenariofile.name)

        args = argparse.Namespace()        
        args.db = dbfile.name
        args.description = None
        args.skip_fragmentation = False
        args.max_broken_bonds = 3
        args.max_water_losses = 1
        args.ms_intensity_cutoff = 0
        args.msms_intensity_cutoff = 0
        args.use_all_peaks = True
        args.adducts = None
        args.max_charge = 1
        args.log = 'debug'
        args.call_back_url = None
        args.scans= 'all'
        args.structure_database = ''
        args.db_options=''
        args.molids = None
        args.ncpus = 1
        args.fast = True
        args.time_limit = None

        self.mc.annotate(args)

        ms = magma.MagmaSession(dbfile.name)
        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset(
                          {
                          'ionisation_mode': 1,
                          'abs_peak_cutoff': 0,
                          'max_ms_level': 5,
                          'mz_precision': 10,
                          'precursor_mz_precision': 0.005,
                          'max_broken_bonds': 3,
                          'ms_intensity_cutoff': 0
                           }, rundata.__dict__)
        scandata = ms.db_session.query(Scan).count()
        self.assertEqual(scandata,5)
        peakdata = ms.db_session.query(Peak).count()
        self.assertEqual(peakdata,19)
        moleculedata = ms.db_session.query(Molecule).count()
        self.assertEqual(moleculedata,68)
        fragmentdata = ms.db_session.query(Fragment).count()
        self.assertEqual(fragmentdata,122)
        hits = ms.db_session.query(Molecule).filter(Molecule.nhits>0).count()
        self.assertEqual(hits,30)

        os.remove(dbfile.name)
