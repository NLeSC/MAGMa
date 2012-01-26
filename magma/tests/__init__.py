import unittest

class TestPeakType(unittest.TestCase):
    def test_init(self):
        import magma
        peaktype = magma.peaktype(123,456,789)
        self.assertEqual(peaktype.mz, 123)

