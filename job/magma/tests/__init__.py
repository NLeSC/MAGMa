import unittest
import mock
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from rdkit import Chem
from rdkit.Chem import AllChem
import magma
from magma.errors import FileFormatError,DataProcessingError
from magma.models import Base, Molecule, Scan, Peak, Fragment, Run

class TestMagmaSession(unittest.TestCase):
    def test_construct_with_new_db(self):
        ms = magma.MagmaSession(':memory:', u'My description')

        self.assertIsInstance(ms, magma.MagmaSession)
        rundata = ms.db_session.query(Run).one()
        self.assertEqual(rundata.description, u'My description')

    def test_construct_with_existing_db(self):
        # create db with run description
        import tempfile, os
        dbfile = tempfile.NamedTemporaryFile(delete=False)
        engine = create_engine('sqlite:///'+dbfile.name)
        Base.metadata.create_all(engine)
        session = sessionmaker(bind=engine)()
        session.add(Run(description=u'My first description'))
        session.commit()
        dbfile.close()

        ms = magma.MagmaSession(dbfile.name, u'My second description')

        rundata = ms.db_session.query(Run).one()
        # first description is not overwritten
        self.assertEqual(rundata.description, u'My first description')

        os.remove(dbfile.name)

    def test_get_structure_engine_default(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        se = ms.get_structure_engine()

        self.assertIsInstance(se, magma.StructureEngine)

        self.assertFalse(hasattr(se,'pubchem_engine'))


    def test_get_structure_engine_custom(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        args = {
                'pubchem_names': True,
                }
        se = ms.get_structure_engine(**args)

        self.assertIsInstance(se, magma.StructureEngine)

        self.assertIsInstance(se.pubchem_engine, magma.PubChemEngine)

    def test_get_ms_data_engine_default(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        se = ms.get_ms_data_engine()

        self.assertIsInstance(se, magma.MsDataEngine)

        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset({
                                      'ionisation_mode': 1,
                                      'abs_peak_cutoff': 1000,
                                      'max_ms_level': 10,
                                      'mz_precision':5.0,
                                      'precursor_mz_precision':0.005
                                       }, rundata.__dict__)

    def test_get_ms_data_engine_custom(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        args = {
                'ionisation_mode': -1,
                'abs_peak_cutoff': 10001,
                'max_ms_level': 101,
                'mz_precision': 0.01,
                'precursor_mz_precision': 0.05
                }
        se = ms.get_ms_data_engine(**args)

        self.assertIsInstance(se, magma.MsDataEngine)

        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset(args, rundata.__dict__)

    def test_get_annotate_engine_default(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        ms.get_ms_data_engine()
        se = ms.get_annotate_engine()

        self.assertIsInstance(se, magma.AnnotateEngine)

        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset({
                 'skip_fragmentation':False,
                 'max_broken_bonds':3,
                 'ms_intensity_cutoff':1e6,
                 'msms_intensity_cutoff':5,
                 'use_all_peaks':False
                                       }, rundata.__dict__)

    def test_get_annotate_engine_custom(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        ms.get_ms_data_engine()
        args = {
                'skip_fragmentation': True,
                'max_broken_bonds': 5,
                'ms_intensity_cutoff': 1e6,
                'msms_intensity_cutoff': 0.2,
                'use_all_peaks': True
                }
        se = ms.get_annotate_engine(**args)

        self.assertIsInstance(se, magma.AnnotateEngine)

        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset(args, rundata.__dict__)

class TestStructureEngine(unittest.TestCase):
    def setUp(self):
        engine = create_engine('sqlite://')
        Base.metadata.create_all(engine)
        self.db_session = sessionmaker(bind=engine)()
        

    def test_construct_new_db(self):
        se = magma.StructureEngine(self.db_session)

        self.assertIsNone(se.call_back_engine)

    def test_add_structure(self):
        se = magma.StructureEngine(self.db_session, u'test')

        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('CCO'))

        molid = se.add_structure(molblock, 'ethanol', 1.0, False)

        mol = self.db_session.query(Molecule).filter(Molecule.molid==molid).one()
        self.assertDictContainsSubset(
                             {
                              'mol': molblock,
                              'refscore': 1.0,
                              'reactionsequence': {},
                              'inchikey14': u'LFQSCWFLJHTTHZ',
                              'formula': u'C2H6O',
                              'predicted': False,
                              'name': u'ethanol',
                              'mim': 46.0418648147,
                              'logp': -0.0014000000000000123
                             },
                             mol.__dict__
                             )

    def test_add_structure_2nd_replace(self):
        se = magma.StructureEngine(self.db_session, u'test')

        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('CCO'))

        se.add_structure(molblock, 'ethanol', 1.0, 0, 'PARENT', 1)
        molid = se.add_structure(molblock, 'ethanol', 2.0, 0, 'PARENT', 1)

        met = self.db_session.query(Molecule).one()

        self.assertDictContainsSubset(
                             {
                              'molid': molid,
                              'refscore': 2.0
                             },
                             met.__dict__
                             )


    def test_add_structure_2nd_ignore(self):
        se = magma.StructureEngine(self.db_session, u'test')

        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('CCO'))

        molid = se.add_structure(molblock, 'ethanol', 2.0, 0, 'PARENT', 1)
        se.add_structure(molblock, 'ethanol', 1.0, 0, 'PARENT', 1)

        met = self.db_session.query(Molecule).one()
        self.assertDictContainsSubset(
                             {
                              'molid': molid,
                              'refscore': 2.0
                             },
                             met.__dict__
                             )

    def test_metabolize(self):
        """ Treat it as blackbox, do not want to know how reactor is called
        just that metabolisation happened
        """
        se = magma.StructureEngine(self.db_session)
        #molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('Oc1cc(CC2OCCC2)cc(O)c1O'))
        parent_molid = se.read_smiles('Oc1cc(CC2OCCC2)cc(O)c1O')
        se.metabolize(1, 'phase1',endpoints=False)
        self.assertGreater(self.db_session.query(Molecule).filter(Molecule.predicted==True).count(), 0)

    def test_metabolize_unknown_metabolite_type(self):
        se = magma.StructureEngine(self.db_session, u'test')
        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('Oc1cc(CC2OCCC2)cc(O)c1O'))
        parent_molid = se.add_structure(
                                        molblock, '5-(3,4,5)-trihydroxyphenyl-g-valerolactone (F,U)',
                                        1.0, 0, 'PARENT', 1
                                        )
        se.metabolize(parent_molid, u'phase1234')

        # TODO metabolize requires at least one SMIRKS query
        # if none given then raise exception

        self.assertEqual(self.db_session.query(Molecule).count(), 1)

    def test_metabolize_bad_molid(self):
        se = magma.StructureEngine(self.db_session, u'test')
        parent_molid = 1
        se.metabolize(parent_molid, u'phase1')

        self.assertEqual(self.db_session.query(Molecule).count(), 0)

    def test_metabolize_all(self):
        se = magma.StructureEngine(self.db_session)
        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('Oc1cc(CC2OCCC2)cc(O)c1O'))
        parent_molid = se.add_structure(
                                        molblock, '5-(3,4,5)-trihydroxyphenyl-g-valerolactone (F,U)',
                                        1.0, 1
                                        )
        se.metabolize = mock.Mock(return_value=set([1]))

        se.metabolize_all(u'phase1')

        se.metabolize.assert_called_with(parent_molid, u'phase1',False)

class TestMsDataEngine(unittest.TestCase):
    def setUp(self):
        engine = create_engine('sqlite://')
        Base.metadata.create_all(engine)
        self.db_session = sessionmaker(bind=engine)()

    def test_store_manual_subtree_corrupt_endlessloop(self):
        mde = magma.MsDataEngine(self.db_session, 1, 1000, 5, 0.001, 0.005, 3)
        # create corrupt manual tree file
        import tempfile, os
        treefile = tempfile.NamedTemporaryFile(delete=False)
        treefile.write("""320.2: 999 (
    123.12    1
    135.2    3
    139.2    3
    145.18    2
    147.04    1
    149.136    39
    149.92    1
    150.32    2
    150.56    1
    150.88    1
    151.28    1
    153.36    1
    155.04    2
    155.296    2
    159.1    1
    161.28    5
    162.64    1
    163.155    43
    167.145    230
    167.808    4
    168.4    5
    169.06    7
    169.68    1
    175.177    8
    175.52    1
    177.12    4
    179.114    92
    179.509    5
    179.74    1
    179.92    2
    180.168    5
    181.04    1
    181.44    1
    189.28    2
    191.12    1
    195.04    3
    195.28    2
    203.235    23
    203.44    3
    205.2    1
    205.68    1
    206.96    4
    207.2    4
    207.36    1
    208.088    13
    219.232    14
    220.88    1
    221.12    4
    223.28    1
    229.12    5
    257.283    254
    257.92    4
    258.32    2
    258.733    3
    259.32    1
    273.12    3
    275.268    97
    282.88    1
    283.23    4
    291.251    9
    291.52    3
    301.242    254
    302    3
    302.4    2
    319.155    999
    319.8    2
)
        """)
        treefile.close()

        with self.assertRaises(FileFormatError) as cm:
            mde.store_manual_tree(treefile.name,0)

        self.assertEqual(str(cm.exception), 'Corrupt Tree format ...')

        os.remove(treefile.name)

    def test_store_manual_tree_unallowed_element(self):
        mde = magma.MsDataEngine(self.db_session, 1, 1000, 5, 0.001, 0.005, 3)
        # create corrupt manual tree file
        import tempfile, os
        treefile = tempfile.NamedTemporaryFile(delete=False)
        treefile.write('C6H5At: 999 (C6H5: 1000, At: 100)')
        treefile.close()

        with self.assertRaises(FileFormatError) as cm:
            mde.store_manual_tree(treefile.name,-1)

        self.assertEqual(str(cm.exception), 'Element not allowed in formula tree: At')

        os.remove(treefile.name)


class TestAnnotateEngine(unittest.TestCase):
    def setUp(self):
        engine = create_engine('sqlite://')
        Base.metadata.create_all(engine)
        self.db_session = sessionmaker(bind=engine)()

    def test_construct_annotate_engine(self):
        with self.assertRaises(DataProcessingError) as cm:
            ae=magma.AnnotateEngine(self.db_session,False,3,1,0,5,True,adducts='Na,K')

        self.assertEqual(str(cm.exception), 'No MS data parameters read.')

    def test_generate_ions(self):
        mde = magma.MsDataEngine(self.db_session, 1, 1000, 5, 0.001, 0.005, 3)
        ae=magma.AnnotateEngine(self.db_session,False,3,1,0,5,True,adducts='Na,K')
        self.assertIsInstance(ae, magma.AnnotateEngine)
        self.assertEqual(ae.ions,[{0: '[M]+'}, {1.0078250321: '[M+H]+', 22.9897692809: '[M+Na]+', 38.96370668: '[M+K]+'}])


class TestMetabolizeEngine(unittest.TestCase):
    def test_metabolize_lumiracoxib_phase1and2(self):
        me = magma.MetabolizeEngine()
        mol=Chem.MolFromSmiles('Clc1cccc(F)c1Nc2ccc(C)cc2CC(=O)O')
        AllChem.Compute2DCoords(mol)
        metabolites=me.metabolize(mol,'phase1',False)
        self.assertEqual(len(metabolites),9)

    def test_glycosidase_endpoints(self):
        molblock="""10078714
  -OEChem-06121220462D

 52 57  0     1  0  0  0  0  0999 V2000
    7.1962   -1.7500    0.0000 O   0  0  0
    8.0622   -4.2500    0.0000 O   0  0  0
    9.7942   -4.2500    0.0000 O   0  0  0
    6.3301   -4.2500    0.0000 O   0  0  0
    6.3301   -0.2500    0.0000 O   0  0  0
    4.5981   -3.2500    0.0000 O   0  0  0
    4.5981   -1.2500    0.0000 O   0  0  0
    8.0622   -6.2500    0.0000 O   0  0  0
    9.7942   -7.2500    0.0000 O   0  0  0
   11.5263   -6.2500    0.0000 O   0  0  0
    2.8660    4.7500    0.0000 O   0  0  0
    3.7320    3.2500    0.0000 O   0  0  0
    5.4641    6.2500    0.0000 O   0  0  0
    3.7320    7.2500    0.0000 O   0  0  0
    5.4641    4.2500    0.0000 O   0  0  0
    8.0622    1.7500    0.0000 O   0  0  0
    8.0622   -1.2500    0.0000 O   0  0  0
    9.8106   -1.2846    0.0000 O   0  0  0
   11.5923    1.7742    0.0000 O   0  0  0
    7.1962   -2.7500    0.0000 C   0  0  1
    6.3301   -3.2500    0.0000 C   0  0  2
    5.4641   -2.7500    0.0000 C   0  0  1
    5.4641   -1.7500    0.0000 C   0  0  2
    8.9282   -5.7500    0.0000 C   0  0  1
    8.9282   -4.7500    0.0000 C   0  0  2
    9.7942   -6.2500    0.0000 C   0  0  1
   10.6603   -5.7500    0.0000 C   0  0  2
    6.3301   -1.2500    0.0000 C   0  0  1
    8.0622   -3.2500    0.0000 C   0  0  0
   10.6603   -4.7500    0.0000 C   0  0  2
    4.5981    5.7500    0.0000 C   0  0  2
    3.7320    6.2500    0.0000 C   0  0  2
    4.5981    4.7500    0.0000 C   0  0  2
    2.8660    5.7500    0.0000 C   0  0  2
   11.5263   -4.2500    0.0000 C   0  0  0
    3.7320    4.2500    0.0000 C   0  0  1
    2.0000    6.2500    0.0000 C   0  0  0
    7.1962    0.2500    0.0000 C   0  0  0
    7.1962    1.2500    0.0000 C   0  0  0
    8.0622   -0.2500    0.0000 C   0  0  0
    4.5981    2.7500    0.0000 C   0  0  0
    6.3301    1.7500    0.0000 C   0  0  0
    8.9282    0.2500    0.0000 C   0  0  0
    4.5981    1.7500    0.0000 C   0  0  0
    5.4641    3.2500    0.0000 C   0  0  0
    8.9282    1.2500    0.0000 C   0  0  0
    5.4641    1.2500    0.0000 C   0  0  0
    6.3301    2.7500    0.0000 C   0  0  0
    9.8222   -0.2847    0.0000 C   0  0  0
    9.8222    1.7847    0.0000 C   0  0  0
   10.7282    0.2292    0.0000 C   0  0  0
   10.7282    1.2708    0.0000 C   0  0  0
  1 20  1  0
  1 28  1  0
 25  2  1  0
  2 29  1  0
  3 25  1  0
  3 30  1  0
 21  4  1  0
 28  5  1  0
  5 38  1  0
 22  6  1  0
 23  7  1  0
 24  8  1  0
 26  9  1  0
 27 10  1  0
 11 34  1  0
 11 36  1  0
 36 12  1  0
 12 41  1  0
 31 13  1  0
 32 14  1  0
 33 15  1  0
 16 39  1  0
 16 46  1  0
 17 40  2  0
 18 49  1  0
 19 52  1  0
 20 21  1  0
 20 29  1  0
 21 22  1  0
 22 23  1  0
 23 28  1  0
 24 25  1  0
 24 26  1  0
 26 27  1  0
 27 30  1  0
 30 35  1  0
 31 32  1  0
 31 33  1  0
 32 34  1  0
 33 36  1  0
 34 37  1  0
 38 39  2  0
 38 40  1  0
 39 42  1  0
 40 43  1  0
 41 44  2  0
 41 45  1  0
 42 47  2  0
 42 48  1  0
 43 46  1  0
 43 49  2  0
 44 47  1  0
 45 48  2  0
 46 50  2  0
 49 51  1  0
 50 52  1  0
 51 52  2  0
M  END
"""
        me = magma.MetabolizeEngine()
        mol=Chem.MolFromMolBlock(molblock)
        metabolites=me.metabolize(mol,'glycosidase',True)
        self.assertEqual(len(metabolites),3)

