import unittest
import mock
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from rdkit import Chem
import magma
from magma.models import Base, Metabolite, Scan, Peak, Fragment, Run

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

        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset({
                                       'metabolism_types': u'phase1,phase2',
                                       'n_reaction_steps': 2
                                       },rundata.__dict__)


    def test_get_structure_engine_custom(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        args = {
                'metabolism_types': u'phase2',
                'n_reaction_steps': 3
                }
        se = ms.get_structure_engine(**args)

        self.assertIsInstance(se, magma.StructureEngine)

        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset(args, rundata.__dict__)

    def test_get_ms_data_engine_default(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        se = ms.get_ms_data_engine()

        self.assertIsInstance(se, magma.MsDataEngine)

        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset({
                                      'abs_peak_cutoff': 1000,
                                      'max_ms_level': 10
                                       }, rundata.__dict__)

    def test_get_ms_data_engine_custom(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        args = {
                'abs_peak_cutoff': 10001,
                'max_ms_level': 101
                }
        se = ms.get_ms_data_engine(**args)

        self.assertIsInstance(se, magma.MsDataEngine)

        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset(args, rundata.__dict__)

    def test_get_annotate_engine_default(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        se = ms.get_annotate_engine()

        self.assertIsInstance(se, magma.AnnotateEngine)

        rundata = ms.db_session.query(Run).one()
        self.assertDictContainsSubset({
                 'ionisation_mode': 1,
                 'skip_fragmentation':False,
                 'max_broken_bonds':4,
                 'ms_intensity_cutoff':1e6,
                 'msms_intensity_cutoff':0.1,
                 'mz_precision':0.001,
                 'precursor_mz_precision':0.005,
                 'use_all_peaks':False
                                       }, rundata.__dict__)

    def test_get_annotate_engine_custom(self):
        ms = magma.MagmaSession(':memory:', u'My description')
        args = {
                'ionisation_mode': -1,
                'skip_fragmentation': True,
                'max_broken_bonds': 5,
                'ms_intensity_cutoff': 1e6,
                'msms_intensity_cutoff': 0.2,
                'mz_precision': 0.01,
                'precursor_mz_precision': 0.05,
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
        se = magma.StructureEngine(self.db_session, u'phase1,phase2', 2)

        self.assertEquals(se.metabolism_types, [ 'phase1', 'phase2'])
        self.assertEquals(se.n_reaction_steps, 2)
        rundata = self.db_session.query(Run).one()
        self.assertDictContainsSubset({
                                       'metabolism_types': u'phase1,phase2',
                                       'n_reaction_steps': 2
                                       }, rundata.__dict__)

    def test_construct_with_existing_db(self):
        self.db_session.add(Run(
                               metabolism_types=u'phase1,phase2',
                               n_reaction_steps=2
                                ))
        self.db_session.commit()

        se = magma.StructureEngine(self.db_session, u'phase2', 3)

        self.assertEquals(se.metabolism_types, [ 'phase1', 'phase2'])
        self.assertEquals(se.n_reaction_steps, 2)
        rundata = self.db_session.query(Run).one()
        self.assertDictContainsSubset({
                                       'metabolism_types': u'phase1,phase2',
                                       'n_reaction_steps': 2
                                       }, rundata.__dict__)

    def test_add_structure(self):
        se = magma.StructureEngine(self.db_session, u'phase1,phase2', 2)

        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('CCO'))

        metid = se.add_structure(molblock, 'ethanol', 1.0, 0, 'PARENT', 1)

        met = self.db_session.query(Metabolite).filter(Metabolite.metid==metid).one()
        self.assertDictContainsSubset(
                             {
                              'mol': molblock,
                              'level': 0,
                              'probability': 1.0,
                              'reactionsequence': u'PARENT',
                              'smiles': u'CCO',
                              'molformula': u'C2H6O',
                              'isquery': True,
                              'origin': u'ethanol',
                              'mim': 46.0418648147,
                              'logp': -0.0014000000000000123
                             },
                             met.__dict__
                             )

    def test_add_structure_2nd_replace(self):
        se = magma.StructureEngine(self.db_session, u'phase1,phase2', 2)

        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('CCO'))

        se.add_structure(molblock, 'ethanol', 1.0, 0, 'PARENT', 1)
        metid = se.add_structure(molblock, 'ethanol', 2.0, 0, 'PARENT', 1)

        met = self.db_session.query(Metabolite).one()

        self.assertDictContainsSubset(
                             {
                              'metid': metid,
                              'probability': 2.0
                             },
                             met.__dict__
                             )


    def test_add_structure_2nd_ignore(self):
        se = magma.StructureEngine(self.db_session, u'phase1,phase2', 2)

        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('CCO'))

        metid = se.add_structure(molblock, 'ethanol', 2.0, 0, 'PARENT', 1)
        se.add_structure(molblock, 'ethanol', 1.0, 0, 'PARENT', 1)

        met = self.db_session.query(Metabolite).one()
        self.assertDictContainsSubset(
                             {
                              'metid': metid,
                              'probability': 2.0
                             },
                             met.__dict__
                             )

    def test_metabolize(self):
        """ Treat it as blackbox, do not want to know how reactor is called
        just that metabolisation happened
        """
        se = magma.StructureEngine(self.db_session, u'phase1', 1)
        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('Oc1cc(CC2OCCC2)cc(O)c1O'))
        parent_metid = se.add_structure(
                                        molblock, '5-(3,4,5)-trihydroxyphenyl-g-valerolactone (F,U)',
                                        1.0, 0, 'PARENT', 1
                                        )

        se.metabolize(parent_metid, u'phase1', 1)

        self.assertGreater(self.db_session.query(Metabolite).filter(Metabolite.level==1).count(), 0)

    def test_metabolize_unknown_metabolite_type(self):
        se = magma.StructureEngine(self.db_session, u'phase1234', 1)
        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('Oc1cc(CC2OCCC2)cc(O)c1O'))
        parent_metid = se.add_structure(
                                        molblock, '5-(3,4,5)-trihydroxyphenyl-g-valerolactone (F,U)',
                                        1.0, 0, 'PARENT', 1
                                        )
        se.metabolize(parent_metid, u'phase1234', 1)

        # TODO reactor requires at least one SMIRKS query
        # if none given then raise exception

        self.assertEqual(self.db_session.query(Metabolite).count(), 1)

    def test_metabolize_bad_metid(self):
        se = magma.StructureEngine(self.db_session, u'phase1', 1)
        parent_metid = 1
        se.metabolize(parent_metid, u'phase1', 1)

        self.assertEqual(self.db_session.query(Metabolite).count(), 0)

    def test_metabolize_all(self):
        se = magma.StructureEngine(self.db_session, u'phase1', 1)
        molblock = Chem.MolToMolBlock(Chem.MolFromSmiles('Oc1cc(CC2OCCC2)cc(O)c1O'))
        parent_metid = se.add_structure(
                                        molblock, '5-(3,4,5)-trihydroxyphenyl-g-valerolactone (F,U)',
                                        1.0, 0, 'PARENT', 1
                                        )
        se.metabolize = mock.Mock()

        se.metabolize_all(u'phase1', 1)

        se.metabolize.assert_called_with(parent_metid, u'phase1', 1)

class TestMsDataEngine(unittest.TestCase):
    def setUp(self):
        engine = create_engine('sqlite://')
        Base.metadata.create_all(engine)
        self.db_session = sessionmaker(bind=engine)()

    def test_store_manual_subtree_corrupt_endlessloop(self):
        mde = magma.MsDataEngine(self.db_session, 1000, 1, 3)
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

        with self.assertRaises(SystemExit) as cm:
            mde.store_manual_tree(treefile.name)

        self.assertEqual(cm.exception.code, 'Corrupt Tree format ...')

        os.remove(treefile.name)


class TestAnnotateEngine(unittest.TestCase):
    pass



