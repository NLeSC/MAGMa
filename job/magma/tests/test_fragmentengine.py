import unittest
import magma.fragmentation_py
import magma.fragmentation_cy

haloperidol = """502
  Mrv0541 09041213592D          

 26 28  0  0  0  0            999 V2000
    5.4436    5.1416    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    2.3645   -5.1416    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    3.3810    2.9981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6500   -1.4290    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7935    0.6335    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.7935    2.2836    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0791    1.8711    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5080    1.8711    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0791    1.0460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5080    1.0460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7935   -0.1915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2060    2.9981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0791   -0.6040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7935    3.7126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0310    2.9981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0791   -1.4290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2060    4.4270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4436    3.7126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3645   -1.8415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0310    4.4270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3645   -2.6665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0791   -3.0791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6500   -3.0791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0791   -3.9041    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6500   -3.9041    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3645   -4.3166    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1 20  1  0  0  0  0
  2 26  1  0  0  0  0
  3  6  1  0  0  0  0
  4 19  2  0  0  0  0
  5  9  1  0  0  0  0
  5 10  1  0  0  0  0
  5 11  1  0  0  0  0
  6  7  1  0  0  0  0
  6  8  1  0  0  0  0
  6 12  1  0  0  0  0
  7  9  1  0  0  0  0
  8 10  1  0  0  0  0
 11 13  1  0  0  0  0
 12 14  2  0  0  0  0
 12 15  1  0  0  0  0
 13 16  1  0  0  0  0
 14 17  1  0  0  0  0
 15 18  2  0  0  0  0
 16 19  1  0  0  0  0
 17 20  2  0  0  0  0
 18 20  1  0  0  0  0
 19 21  1  0  0  0  0
 21 22  2  0  0  0  0
 21 23  1  0  0  0  0
 22 24  1  0  0  0  0
 23 25  2  0  0  0  0
 24 26  2  0  0  0  0
 25 26  1  0  0  0  0
M  END
"""


class TestFragmentEnginePython(unittest.TestCase):
    """Run test with python version of FragmentEngine"""
    def setUp(self):
        self.FragmentEngine = magma.fragmentation_py.FragmentEngine

    def test_it(self):
        fe = self.FragmentEngine(mol=haloperidol,
                            max_broken_bonds=3,
                            max_water_losses=1,
                            ionisation_mode=1,
                            skip_fragmentation=0,
                            molcharge=0
                            )
        self.assertEqual(fe.get_natoms(),26)
        self.assertEqual(fe.accepted(),True)
        nfrags = fe.generate_fragments()
        self.assertEqual(nfrags,768)
        atomstring,atomlist,formula,inchikey=fe.get_fragment_info(1,0) # fragment is atom 1
        self.assertEqual(formula, 'Cl')
        atomstring,atomlist,formula,inchikey=fe.get_fragment_info(36,0) # fragment is atom 3 and 6
        self.assertEqual(formula, 'CHO')
        self.assertEqual(atomlist, [2,5])
        self.assertEqual(atomstring, '2,5')
        fragments=fe.find_fragments(181.0, 2^26-1,1.0,0.1)
        self.assertEqual(len(fragments[0]),5)


class TestFragmentEngineCython(TestFragmentEnginePython):
    """Run test with cython version of FragmentEngine"""
    def setUp(self):
        self.FragmentEngine = magma.fragmentation_cy.FragmentEngine
