import unittest

import numpy as np

try:
    from . import spheroid_performance as sph
except ImportError:
    import spheroid_performance as sph


class TestSpheroidPerformanceExample(unittest.TestCase):

    def test_centered_sphere_scales_to_unit_box_guard_band(self):
        case = sph.CASES["centered_sphere"]

        np.testing.assert_allclose(case.semi_axes, (0.49, 0.49, 0.49))
        self.assertAlmostEqual(case.scale, 0.49)
        np.testing.assert_allclose(case.extents, ((0.01, 0.99), (0.01, 0.99), (0.01, 0.99)))
        np.testing.assert_allclose(case.xy_slice_axes(0.5), (0.49, 0.49))
        np.testing.assert_allclose(case.xy_slice_axes(0.01), (0.0, 0.0))

    def test_centered_sphere_slice_and_track_are_hand_calculable(self):
        case = sph.CASES["centered_sphere"]
        model = case.model()

        product = model.sample_slice_xy(case.center[2], 5, 5)
        np.testing.assert_array_equal(
            product.region_id,
            [
                [1, 0, 0, 0, 1],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [1, 0, 0, 0, 1],
            ],
        )

        tracked = model.track(
            vectors=sph.DIRECTION_SETS["x"], spacing=10.0, minimum_tracks=1
        )

        np.testing.assert_array_equal(tracked.segments.region_id, [1, 0, 1])
        np.testing.assert_allclose(tracked.segments.length, [0.01, 0.98, 0.01])
        np.testing.assert_allclose(case.centerline_lengths(0), [0.01, 0.98, 0.01])

    def test_offcenter_cases_keep_each_coordinate_axis_visible(self):
        cases = [
            ("offcenter_x", "x", 0, [0.01, 0.68, 0.31]),
            ("offcenter_y", "y", 1, [0.31, 0.68, 0.01]),
            ("offcenter_z", "z", 2, [0.01, 0.68, 0.31]),
        ]

        for case_name, direction_name, axis, lengths in cases:
            with self.subTest(case=case_name):
                case = sph.CASES[case_name]
                tracked = case.model().track(
                    vectors=sph.DIRECTION_SETS[direction_name],
                    spacing=10.0,
                    minimum_tracks=1,
                )

                np.testing.assert_array_equal(tracked.segments.region_id, [1, 0, 1])
                np.testing.assert_allclose(tracked.segments.length, lengths)
                np.testing.assert_allclose(case.centerline_lengths(axis), lengths)

    def test_benchmark_case_reports_timing_counts_and_volume_error(self):
        result = sph.benchmark_case(
            sph.CASES["oblong_x"],
            directions=sph.DIRECTION_SETS["xyz"],
            spacing=10.0,
            slice_cells=8,
        )

        self.assertEqual(result["name"], "oblong_x")
        self.assertEqual(result["slice_cells"], 8)
        self.assertGreater(result["tracks"], 0)
        self.assertGreater(result["segments"], 0)
        self.assertGreater(result["timings"]["total_s"], 0.0)
        self.assertIn(sph.SPHEROID_REGION, result["tracked_region_volumes"])
        self.assertAlmostEqual(result["exact_total_volume"], 1.0)
        self.assertIn("total_volume_error", result)
        np.testing.assert_allclose(
            result["analytic_slice_xy_bounds"], [0.01, 0.99, 0.255, 0.745]
        )
        np.testing.assert_allclose(
            result["analytic_projection_xy_bounds"], [0.01, 0.99, 0.255, 0.745]
        )
        self.assertAlmostEqual(
            result["exact_spheroid_volume"],
            4.0 * np.pi * 0.49 * 0.245 * 0.1225 / 3.0,
        )


if __name__ == "__main__":
    unittest.main()
