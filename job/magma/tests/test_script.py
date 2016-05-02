import unittest
import argparse
import magma.script
import tempfile, os
import pkg_resources
from StringIO import StringIO
from magma.models import Base, Molecule, Scan, Peak, Fragment, Run

class TestMagmaCommand(unittest.TestCase):

    def setUp(self):
        self.mc = magma.script.MagmaCommand()

    def test_version(self):

        self.assertEqual(self.mc.version(), '1.0')

    def test_theogallin_example(self):
        mzxmlfile = pkg_resources.resource_filename('magma', "tests/theogallin.mzXML")
        sdfile = pkg_resources.resource_filename('magma', "tests/theogallin.sdf")
        dbfile = tempfile.NamedTemporaryFile(delete=False)

        args = argparse.Namespace()
        args.db = dbfile.name
        args.ms_data = mzxmlfile
        args.description = 'Theogallin_example'
        args.ionisation_mode = -1
        args.abs_peak_cutoff = 10000
        args.mz_precision = 5
        args.mz_precision_abs = 0.001
        args.precursor_mz_precision = 0.005
        args.max_ms_level = 5
        args.ms_data_format = 'mzxml'
        args.scan=None
        #args.ms_data = argparse.Namespace(name='bogus.mzxml')
        args.log = 'debug'
        args.call_back_url = None
        args.time_limit = None

        self.mc.read_ms_data(args)

        args = argparse.Namespace()        
        args.db = dbfile.name
        args.description = None
        args.log = 'debug'
        args.pubchem_names = False
        args.structure_format = 'sdf'
        args.mass_filter = 9999
        args.structures = sdfile
        self.mc.add_structures(args)

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
        args.fast = False
        args.time_limit = None

        self.mc.annotate(args)

        ms = magma.MagmaSession(dbfile.name)
        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset(
                          {
                          'ionisation_mode': -1,
                          'abs_peak_cutoff': 10000,
                          'max_ms_level': 5,
                          'mz_precision': 5.0,
                          'precursor_mz_precision': 0.005,
                          'max_broken_bonds': 3,
                          'ms_intensity_cutoff': 0
                           }, rundata.__dict__)
        scandata = ms.db_session.query(Scan).count()
        self.assertEqual(scandata,4)
        peakdata = ms.db_session.query(Peak).count()
        self.assertEqual(peakdata,32)
        moleculedata = ms.db_session.query(Molecule).count()
        self.assertEqual(moleculedata,1)
        fragmentdata = ms.db_session.query(Fragment).count()
        self.assertEqual(fragmentdata,13)

        os.remove(dbfile.name)

    def test_chlorogenic_acid_example_without_fast_option(self):
        treefile = tempfile.NamedTemporaryFile(delete=False)
        dbfile = tempfile.NamedTemporaryFile(delete=False)

        args = argparse.Namespace()
        args.db = dbfile.name
        args.ms_data = treefile.name
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

        # allow,but not require, commas at the end of a line
        treefile.write("""353.087494: 69989984 (
    191.055756: 54674544 (
        85.029587: 2596121,
        93.034615: 1720164
        109.029442: 917026
        111.045067: 1104891 (
\t          81.034691: 28070,
            83.014069:7618,
            83.050339:25471,
            93.034599: 36300,
            96.021790: 8453
            )
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
        args.scans= '1'
        args.structure_database = 'hmdb'
        args.db_options=pkg_resources.resource_filename('magma', "tests/HMDB_MAGMa_test.db")
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
        self.assertEqual(moleculedata,5)
        fragmentdata = ms.db_session.query(Fragment).count()
        self.assertEqual(fragmentdata,87)

        sdfile = tempfile.NamedTemporaryFile(delete=False)
        args = argparse.Namespace()        
        args.db = dbfile.name
        args.description = None
        args.filename = sdfile.name
        args.assigned = None
        args.output_format = 'sdf'
        self.mc.export_structures(args)
        
        scores=''
        sl=False
        for line in open(sdfile.name,'r'):
            if sl:
                scores+=line
                sl=False
            elif line == '> <score>\n':
                sl=True
        self.assertEqual(scores,'1.77587\n1.78058\n1.78367\n7.05795\n7.05795\n')

        os.remove(sdfile.name)
        os.remove(dbfile.name)

    def test_light_glutathion_mgf(self):
        treefile = tempfile.NamedTemporaryFile(delete=False)

        args = argparse.Namespace()
        args.ms_data = treefile.name
        args.description = 'Example'
        args.ionisation_mode = 1
        args.abs_peak_cutoff = 0
        args.mz_precision = 5
        args.mz_precision_abs = 0.001
        args.precursor_mz_precision = 0.005
        args.max_ms_level = 5
        args.ms_data_format = 'mgf'
        args.log = 'debug'
        args.call_back_url = None

        # allow,but not require, commas at the end of a line
        treefile.write("""BEGIN IONS
TITLE=CASMI 2014, Challenge 9
PEPMASS=308.0912 100.0
16.0165 3.2
144.0114 6.3
162.0219 40.2
179.0485 100.0
233.0590 21.6
290.0802 5.1
END IONS
""")
        treefile.close()

        args.skip_fragmentation = False
        args.max_broken_bonds = 3
        args.max_water_losses = 1
        args.adducts = None
        args.max_charge = 1
        args.ncpus = 1
        args.slow = False
        args.time_limit = None
        args.structure_database = ""

        args.structure_database = 'hmdb'
        args.db_options=pkg_resources.resource_filename('magma', "tests/HMDB_MAGMa_test.db")
        args.read_molecules = 'NC(CCC(=O)NC(CS)C(=O)NCC(=O)O)C(=O)O'
        args.output_format = 'smiles'
        out = StringIO()
        self.mc.light(args, out)
        
        self.assertEqual(out.getvalue(),u'NC(CCC(=O)NC(CS)C(=O)NCC(=O)O)C(=O)O score=1.22932 name= refscore=None formula=C10H17N3O6S mim=307.083805984\n')

        os.remove(treefile.name)
        

    def test_JWH015_example_with_fast_option(self):
        dbfile = tempfile.NamedTemporaryFile(delete=False)
        treefile = tempfile.NamedTemporaryFile(delete=False)
        args = argparse.Namespace()
        args.db = dbfile.name
        args.ms_data = treefile.name
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
        scenariofile.write('phase1_selected,1\nmass_filter,500\nphase2_selected,2')
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
        self.assertEqual(moleculedata,39)
        fragmentdata = ms.db_session.query(Fragment).count()
        self.assertEqual(fragmentdata,85)
        hits = ms.db_session.query(Molecule).filter(Molecule.nhits>0).count()
        self.assertEqual(hits,21)

        sdfile = tempfile.NamedTemporaryFile(delete=False)
        args = argparse.Namespace()        
        args.db = dbfile.name
        args.description = None
        args.filename = sdfile.name
        args.assigned = None
        args.output_format = 'sdf'
        self.mc.export_structures(args)
        nonmatched=0
        hl=False
        for line in open(sdfile.name,'r'):
            if hl:
                if line=='0\n':
                    nonmatched+=1
                hl=False
            elif line == '> <nhits>\n':
                hl=True
        self.assertEqual(nonmatched,18)
        os.remove(sdfile.name)

        #assign parent to first peak in data
        peak=ms.db_session.query(Peak).first()
        peak.assigned_molid=1
        ms.db_session.add(peak)
        ms.db_session.commit()

        sdfile = tempfile.NamedTemporaryFile(delete=False)
        args = argparse.Namespace()        
        args.db = dbfile.name
        args.description = None
        args.filename = sdfile.name
        args.assigned = True
        self.mc.export_structures(args)
        l=''
        f=open(sdfile.name,'r')
        while l != '> <mz>\n':
            l=f.readline()
        l=f.readline()
        self.assertEqual(l, '328.1696\n')
        os.remove(sdfile.name)
        
        os.remove(dbfile.name)
        
