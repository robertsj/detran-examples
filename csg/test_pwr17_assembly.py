import tempfile
import unittest
from pathlib import Path

import numpy as np

try:
    from . import pwr17_assembly as pwr
except ImportError:
    import pwr17_assembly as pwr


class TestPWR17AssemblyCSG(unittest.TestCase):

    def test_example_pixel_map_defaults_are_five_times_the_scaffold_resolution(self):
        self.assertEqual(pwr.SLICE_RESOLUTION_SCALE, 5)
        self.assertEqual(
            pwr.DEFAULT_HOMOGENEOUS_SLICE_CELLS,
            5 * pwr.BASE_HOMOGENEOUS_SLICE_CELLS,
        )
        self.assertEqual(
            pwr.DEFAULT_PIN_SLICE_CELLS,
            5 * pwr.BASE_PIN_SLICE_CELLS,
        )
        self.assertEqual(
            pwr.DEFAULT_ASSEMBLY_SLICE_CELLS,
            5 * pwr.BASE_ASSEMBLY_SLICE_CELLS,
        )

    def test_homogeneous_exit_ramp_is_a_single_region_box(self):
        model = pwr.homogeneous_model()

        self.assertEqual(model.region_at(1.0, 1.0).name, "moderator")
        self.assertEqual(model.region_at(model.width_x + 0.1, 1.0).region_id, pwr.OUTSIDE)

        product = model.sample_slice(4, 3)

        self.assertEqual(product.region_id.shape, (3, 4))
        np.testing.assert_array_equal(product.region_id, np.zeros((3, 4), dtype=int))

    def test_fuel_pin_2d_walks_from_fuel_to_gap_to_clad_to_moderator(self):
        model = pwr.fuel_pin_2d_model()
        c = 0.5 * model.spec.pin_pitch

        self.assertEqual(model.region_at(c, c).name, "uo2_fuel")
        self.assertEqual(
            model.region_at(c + 0.414, c).name,
            "helium_gap",
        )
        self.assertEqual(
            model.region_at(c + 0.45, c).name,
            "zircaloy_clad",
        )
        self.assertEqual(
            model.region_at(c + 0.60, c).name,
            "moderator",
        )

    def test_assembly_slice_has_standard_24_guides_and_center_instrument(self):
        model = pwr.assembly_2d_model(control_position=0.0)

        self.assertEqual(len(pwr.GUIDE_TUBE_LOCATIONS), 24)
        self.assertEqual(model.pin_kind(8, 8), "instrument")
        self.assertEqual(model.pin_kind(5, 2), "guide")
        self.assertEqual(model.pin_kind(0, 0), "fuel")

        x, y = model.spec.pin_center(5, 2)
        self.assertEqual(model.region_at(x, y).name, "guide_tube_water")
        x, y = model.spec.pin_center(8, 8)
        self.assertEqual(model.region_at(x, y).name, "instrument_void")

    def test_control_position_changes_guide_tube_material_with_height(self):
        inserted = pwr.assembly_3d_model(control_position=1.0)
        withdrawn = pwr.assembly_3d_model(control_position=0.0)
        x, y = inserted.spec.pin_center(5, 2)

        self.assertEqual(inserted.region_at(x, y, z=50.0).name, "control_aic")
        self.assertEqual(inserted.region_at(x, y, z=200.0).name, "control_b4c")
        self.assertEqual(withdrawn.region_at(x, y, z=200.0).name, "guide_tube_water")

        half_inserted = pwr.assembly_3d_model(control_position=0.5)
        bottom, top = half_inserted.control_poison_bounds()

        self.assertGreater(bottom, inserted.control_poison_bounds()[0])
        self.assertLess(top, withdrawn.control_poison_bounds()[1])

    def test_spacer_grids_are_optional_axial_detail(self):
        with_grids = pwr.assembly_3d_model(include_spacer_grids=True)
        no_grids = pwr.assembly_3d_model(include_spacer_grids=False)
        ix, iy = 0, 0
        cx, cy = with_grids.spec.pin_center(ix, iy)
        x = cx + 0.60

        self.assertEqual(with_grids.region_at(x, cy, z=75.2).name, "zircaloy_spacer_grid")
        self.assertEqual(no_grids.region_at(x, cy, z=75.2).name, "moderator")
        self.assertEqual(with_grids.region_at(x, cy, z=90.0).name, "moderator")

    def test_slice_product_can_be_rendered_as_region_colored_svg(self):
        model = pwr.assembly_3d_model(control_position=1.0)
        product = model.sample_slice(85, 85, z=200.0)

        self.assertEqual(product.region_id.shape, (85, 85))
        self.assertIn(pwr.FUEL, product.region_id)
        self.assertIn(pwr.CONTROL_B4C, product.region_id)

        try:
            svg = pwr.slice_to_svg(product, width=240)
        except ImportError:
            self.skipTest("pydetran visualization helpers are not importable")
        self.assertIn('data-detran-product="slice"', svg)
        self.assertIn('data-region-id="9"', svg)

    def test_slice_product_can_be_written_as_a_pixel_map_png(self):
        try:
            from PIL import Image
        except ImportError:
            self.skipTest("Pillow is optional and not installed")

        model = pwr.assembly_3d_model(control_position=1.0)
        product = model.sample_slice(17, 19, z=200.0)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "assembly.png"
            pwr.slice_to_png(product, path)
            with Image.open(path) as image:
                self.assertEqual(image.size, (17, 19))
                self.assertEqual(image.mode, "RGB")

    def test_ray_trace_product_collapses_path_into_region_segments(self):
        model = pwr.fuel_pin_2d_model()

        rays = model.trace_rays(angle_degrees=0.0, number_rays=1, sample_count=600)

        self.assertEqual(rays.number_rays, 1)
        self.assertGreater(rays.number_segments, 4)
        self.assertEqual(rays.segment_track_index.tolist(), [0] * rays.number_segments)
        self.assertEqual(rays.segment_local_index.tolist(), list(range(rays.number_segments)))
        self.assertAlmostEqual(rays.segment_s_in[0], 0.0)
        self.assertAlmostEqual(rays.segment_s_out[-1], rays.track_length[0])
        self.assertIn(pwr.FUEL, rays.segment_region_id)
        self.assertIn(pwr.GAP, rays.segment_region_id)
        self.assertIn(pwr.CLAD, rays.segment_region_id)
        self.assertIn(pwr.MODERATOR, rays.segment_region_id)

    def test_ray_trace_svg_renders_real_vector_segments(self):
        model = pwr.assembly_3d_model(control_position=1.0)
        rays = model.trace_rays(
            z=200.0,
            angle_degrees=12.0,
            number_rays=9,
            sample_count=300,
        )

        try:
            svg = pwr.ray_trace_to_svg(rays, width=300)
        except ImportError:
            self.skipTest("pydetran visualization helpers are not importable")
        self.assertIn('data-detran-product="tracks"', svg)
        self.assertIn('data-region-id="9"', svg)
        self.assertGreater(svg.count("<line "), 9)

    def test_example_writer_emits_ray_files_for_each_requested_angle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                written = pwr.write_example_outputs(
                    tmpdir,
                    nx=5,
                    ny=5,
                    hdf5=False,
                    number_rays=3,
                    ray_samples=60,
                    ray_angle=(0.0, 45.0),
                )
            except ImportError:
                self.skipTest("pydetran visualization helpers are not importable")

            names = {path.name for path in written}

        self.assertIn("assembly_2d.png", names)
        self.assertNotIn("assembly_2d.svg", names)
        self.assertIn("assembly_3d_midplane_rays_angle_0.svg", names)
        self.assertIn("assembly_3d_spacer_grid_rays_angle_0.svg", names)
        self.assertIn("assembly_3d_midplane_rays_angle_45.svg", names)
        self.assertIn("assembly_3d_spacer_grid_rays_angle_45.svg", names)

    def test_hdf5_export_records_source_and_slice_products(self):
        try:
            import h5py
        except ImportError:
            self.skipTest("h5py is optional and not installed")

        model = pwr.assembly_3d_model(control_position=0.75)
        product = model.sample_slice(10, 12, z=75.2)
        rays = model.trace_rays(z=75.2, number_rays=5, sample_count=100)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "pwr17.h5"
            pwr.write_hdf5(path, model, {"spacer": product}, {"spacer_rays": rays})
            with h5py.File(path, "r") as h5:
                source = h5["geometry/source"]
                self.assertEqual(source.attrs["backend"], "python-csg")
                self.assertAlmostEqual(source.attrs["control_position"], 0.75)
                self.assertIn("regions", h5["geometry"])
                self.assertEqual(
                    h5["geometry/products/slices/spacer/region_id"].shape,
                    (12, 10),
                )
                self.assertIn(
                    "zircaloy_spacer_grid",
                    [
                        name.decode("utf-8") if isinstance(name, bytes) else name
                        for name in h5["geometry/regions/name"][()]
                    ],
                )
                self.assertEqual(
                    h5["geometry/products/rays/spacer_rays/tracks/enter"].shape,
                    (5, 3),
                )
                self.assertIn(
                    pwr.SPACER_ZIRC,
                    h5["geometry/products/rays/spacer_rays/segments/region_id"][()],
                )


if __name__ == "__main__":
    unittest.main()
