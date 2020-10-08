import unittest
import biblis2d, iaea2d, iaea3d, lra2d

class TestDiffusionBenchmarks(unittest.TestCase):

    def test_iaea2d(self):
        got = iaea2d.run()
        expect = 1.02959
        self.assertAlmostEqual(got, expect, 3)

    def test_biblis2d(self):
        got = biblis2d.run()
        expect = 1.02513
        self.assertAlmostEqual(got, expect, 3)

    def test_lra2d(self):
        got = lra2d.run()
        expect = 0.9963
        self.assertAlmostEqual(got, expect, 3)

    def test_iaea3d(self):
        got = iaea3d.run()
        expect = 1.029096
        self.assertAlmostEqual(got, expect, 3)

if __name__ == '__main__':
    unittest.main()
