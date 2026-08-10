import unittest

try:
    from . import run_moc
except ImportError:
    import run_moc


class TestKobayashiMOC(unittest.TestCase):

    def test_straight_duct_runs_native_moc_smoke_cases(self):
        for equation in ("scmoc", "ddmoc"):
            with self.subTest(equation=equation):
                got = run_moc.run(
                    case=2,
                    equation=equation,
                    tolerance=1.0e-3,
                    max_iters=5,
                )

                self.assertEqual(got["case"], 2)
                self.assertEqual(got["equation"], equation)
                self.assertEqual(got["cells"], 360)
                self.assertEqual(got["source_cells"], 1)
                self.assertTrue(got["flux"]["finite"])
                self.assertGreater(got["flux"]["max"], 0.0)
                self.assertGreater(got["sweeps"], 0)


if __name__ == "__main__":
    unittest.main()
