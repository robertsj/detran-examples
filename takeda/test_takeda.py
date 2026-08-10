import unittest

try:
    from . import run_moc, takeda_inserted, takeda_withdrawn
except ImportError:
    import run_moc, takeda_inserted, takeda_withdrawn

class TestTakedaBenchmarks(unittest.TestCase):

    def test_inserted(self):
        got = takeda_inserted.run(write_outputs=False)
        expect = 0.9624
        self.assertAlmostEqual(got, expect, 2)

    def test_withdrawn(self):
        got = takeda_withdrawn.run(write_outputs=False)
        expect = 0.9780
        self.assertAlmostEqual(got, expect, 2)

    def test_native_moc_smoke(self):
        for equation in ("scmoc", "ddmoc"):
            with self.subTest(equation=equation):
                got = run_moc.run(
                    case="inserted",
                    equation=equation,
                    mesh_divisor=10,
                    polar=1,
                    azimuth=1,
                    spacing=10.0,
                    inner_tolerance=1.0e-3,
                    outer_tolerance=1.0e-3,
                    eigen_tolerance=1.0e-3,
                    inner_iters=10,
                    outer_iters=10,
                    eigen_iters=3,
                    monitor_level=0,
                )
                self.assertTrue(got["finite"])
                self.assertGreater(got["keff"], 0.0)
                self.assertEqual(got["cells"], 125)
                self.assertTrue(all(group["finite"] for group in got["flux"]))

if __name__ == '__main__':
    unittest.main()
