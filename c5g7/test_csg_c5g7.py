import inspect
import unittest
from argparse import Namespace

import detran as dt

from .csg_c5g7 import (
  ASSEMBLY_PITCH,
  CORE_PITCH,
  PIN_PITCH,
  PIN_RADIUS,
  PinDiscretization,
  build_csg_overlay_mesh,
  build_assembly,
  build_core,
  build_pin,
  summarize_case,
)


class TestC5G7CSG(unittest.TestCase):

  def test_default_uo2_pin_has_fuel_and_moderator_sectors(self):
    case = build_pin(0)
    self.assertEqual(case.model.number_regions, 16)
    self.assertEqual(len(case.regions), case.model.number_regions)
    self.assertEqual(case.model.material_at(0.5 * PIN_PITCH, 0.5 * PIN_PITCH), 0)
    self.assertEqual(case.model.material_at(0.02, 0.02), 6)
    self.assertEqual(case.model.material_at(0.5 * PIN_PITCH + 0.5 * PIN_RADIUS, 0.5 * PIN_PITCH), 0)

  def test_default_moderator_pin_is_moderator_only(self):
    case = build_pin(6)
    self.assertEqual(case.model.number_regions, 8)
    self.assertEqual({record.material for record in case.regions}, {6})
    self.assertEqual(case.model.material_at(0.5 * PIN_PITCH, 0.5 * PIN_PITCH), 6)
    self.assertEqual(case.model.material_at(PIN_PITCH - 0.02, PIN_PITCH - 0.02), 6)

  def test_radial_refinement_changes_region_count(self):
    case = build_pin(
      0,
      PinDiscretization(
        fuel_radial_divisions=2,
        moderator_radial_divisions=3,
        azimuthal_slices=8,
      ),
    )
    self.assertEqual(case.model.number_regions, 40)

  def test_uo2_assembly_preserves_guide_and_fission_locations(self):
    case = build_assembly(0)
    self.assertAlmostEqual(case.model.width_x, ASSEMBLY_PITCH)
    guide_x = (5 + 0.5) * PIN_PITCH
    guide_y = (2 + 0.5) * PIN_PITCH
    fission_x = (8 + 0.5) * PIN_PITCH
    fission_y = (8 + 0.5) * PIN_PITCH
    self.assertEqual(case.model.material_at(guide_x, guide_y), 4)
    self.assertEqual(case.model.material_at(fission_x, fission_y), 5)
    self.assertEqual(case.model.material_at(0.5 * PIN_PITCH, 0.5 * PIN_PITCH), 0)

  def test_core_dimensions_and_region_metadata_match(self):
    case = build_core()
    self.assertAlmostEqual(case.model.width_x, CORE_PITCH)
    self.assertAlmostEqual(case.model.width_y, CORE_PITCH)
    self.assertEqual(case.model.number_regions, len(case.regions))
    self.assertGreater(case.model.number_regions, 0)

  def test_native_csg_moc_smoke_runs_for_pin(self):
    args = Namespace(
      track=False,
      solve_csg_moc=True,
      csg_overlay_mesh_variant="none",
      csg_overlay_pin_mesh=7,
      solve_mesh_moc=False,
      moc_polar=1,
      moc_azimuth=2,
      moc_quadrature="u-u",
      moc_spacing=2.0,
      moc_equation="scmoc",
      moc_inner_tolerance=1.0e-5,
      moc_outer_tolerance=1.0e-5,
      moc_eigen_tolerance=1.0e-3,
      moc_inner_iters=20,
      moc_outer_iters=20,
      moc_eigen_iters=10,
      moc_monitor_level=0,
      moc_print_interval=10,
    )

    result = summarize_case(build_pin(0), args=args)
    transport = result["transport"]

    self.assertEqual(result["transport_backend"], "native_csg_moc")
    self.assertEqual(transport["cells"], result["number_regions"])
    self.assertGreater(transport["tracks"], 0)
    self.assertGreater(transport["segments"], 0)
    self.assertTrue(transport["finite"])
    self.assertGreater(transport["keff"], 0.0)

  def test_centered_overlay_mesh_metadata_uses_cartesian_planes_only(self):
    mesh, metadata = build_csg_overlay_mesh(
      "pin",
      0,
      variant="centered",
      pin_mesh=2,
    )

    self.assertTrue(metadata["enabled"])
    self.assertEqual(metadata["variant"], "centered")
    self.assertEqual(mesh.number_cells(), 4)
    self.assertIn("materials remain true analytic CSG materials", metadata["material_map_note"])

  def test_native_csg_moc_can_refine_fsrs_with_centered_overlay(self):
    if "overlay_mesh" not in inspect.signature(dt.CSGEigenvalueProblem.__init__).parameters:
      self.skipTest("pydetran CSGEigenvalueProblem overlay_mesh is not available")

    args = Namespace(
      track=False,
      solve_csg_moc=True,
      csg_overlay_mesh_variant="centered",
      csg_overlay_pin_mesh=2,
      solve_mesh_moc=False,
      moc_polar=1,
      moc_azimuth=1,
      moc_quadrature="u-u",
      moc_spacing=2.0,
      moc_equation="scmoc",
      moc_inner_tolerance=1.0e-5,
      moc_outer_tolerance=1.0e-5,
      moc_eigen_tolerance=1.0e-3,
      moc_inner_iters=20,
      moc_outer_iters=20,
      moc_eigen_iters=10,
      moc_monitor_level=0,
      moc_print_interval=10,
    )

    result = summarize_case(build_pin(0), args=args)
    transport = result["transport"]

    self.assertEqual(result["transport_backend"], "native_csg_moc")
    self.assertEqual(transport["fsr_model"], "true_csg_regions_split_by_cartesian_overlay")
    self.assertGreater(transport["cells"], result["number_regions"])
    self.assertEqual(transport["overlay_mesh"]["variant"], "centered")
    self.assertTrue(transport["overlay_product"]["has_region_mesh_maps"])
    self.assertGreater(transport["overlay_product"]["unique_overlay_cells"], 1)
    self.assertTrue(transport["finite"])


if __name__ == "__main__":
  unittest.main()
