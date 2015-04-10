import unittest
import magma.fragmentation_py
import magma.fragmentation_cy

alcohol = """
  Mrv0541 04051115152D

  3  2  0  0  0  0            999 V2000
   -0.6888    0.1371    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0257   -0.2754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7402    0.1371    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
"""


class TestFragmentEnginePython(unittest.TestCase):
    """Run test with python version of FragmentEngine"""
    def setUp(self):
        self.FragmentEngine = magma.fragmentation_py.FragmentEngine

    def test_it(self):
        fe = self.FragmentEngine(mol=alcohol,
                            max_broken_bonds=3,
                            max_water_losses=2,
                            ionisation_mode=1,
                            skip_fragmentation=0,
                            molcharge=0
                            )
        assert fe.get_natoms() == 3


class TestFragmentEngineCython(TestFragmentEnginePython):
    """Run test with cython version of FragmentEngine"""
    def setUp(self):
        self.FragmentEngine = magma.fragmentation_cy.FragmentEngine
