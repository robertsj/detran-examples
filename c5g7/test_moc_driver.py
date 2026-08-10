import unittest

try:
    from .moc_driver import run_moc_case
except ImportError:
    from moc_driver import run_moc_case


class TestC5G7MOCDriver(unittest.TestCase):

    def test_pin_driver_runs_ddmoc_smoke_case(self):
        got = run_moc_case(
            "pin",
            equation="ddmoc",
            index=0,
            pin_mesh=1,
            plot_dir=None,
            spacing=2.0,
            inner_iters=5,
            outer_iters=5,
            eigen_iters=3,
            monitor_level=0,
        )

        self.assertEqual(got["case"], "pin")
        self.assertEqual(got["equation"], "ddmoc")
        self.assertEqual(got["cells"], 1)
        self.assertTrue(got["finite"])
        self.assertGreater(got["keff"], 0.0)
        self.assertGreater(got["sweeps"], 0)


if __name__ == "__main__":
    unittest.main()
