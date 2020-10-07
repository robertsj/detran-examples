import unittest
import takeda_inserted, takeda_withdrawn

class TestDiffusionBenchmarks(unittest.TestCase):

    def test_inserted(self):
        got = takeda_inserted.run()
        expect = 0.9624
        self.assertAlmostEqual(got, expect, 2)

    def test_withdrawn(self):
        got = takeda_withdrawn.run()
        expect = 0.9780
        self.assertAlmostEqual(got, expect, 2)

if __name__ == '__main__':
    unittest.main()
