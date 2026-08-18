import math
import unittest

from .csg_c5g7 import ASSEMBLY_PIN_MAPS, CORE_ASSEMBLY_MAP
from .material_c5g7 import get_density_corrected_materials, get_materials
from .sn_mesh import (
  ASSEMBLY_DIMENSION,
  CORE_DIMENSION,
  PIN_PITCH,
  PIN_RADIUS,
  center_sampled_material_map_for_pin,
  center_sampled_mesh_recipe,
  core_pin_map,
  expanded_materials,
  material_map_for_pin,
  optimized_pin_recipe_7,
  repeated_edges,
)


class TestC5G7SNMesh(unittest.TestCase):

  def test_optimized_7x7_recipe_preserves_pin_area(self):
    recipe = optimized_pin_recipe_7()
    self.assertEqual(recipe.number, 7)
    self.assertEqual(recipe.central_cell_count, 13)
    self.assertAlmostEqual(recipe.edges[0], 0.0)
    self.assertAlmostEqual(recipe.edges[-1], PIN_PITCH)
    self.assertAlmostEqual(sum(recipe.widths), PIN_PITCH)
    self.assertLess(abs(recipe.central_area_error), 1.0e-12)
    self.assertLess(abs(recipe.moderator_area_error), 1.0e-12)
    self.assertAlmostEqual(recipe.central_area_target, math.pi * PIN_RADIUS * PIN_RADIUS)
    self.assertLessEqual(
      recipe.symmetric_difference_fraction,
      recipe.legacy_symmetric_difference_fraction + 1.0e-12,
    )

  def test_optimized_7x7_widths_are_symmetric_and_positive(self):
    recipe = optimized_pin_recipe_7()
    self.assertTrue(all(width > 0.0 for width in recipe.widths))
    self.assertEqual(recipe.widths, tuple(reversed(recipe.widths)))
    self.assertGreater(recipe.optimized_l_over_r, 1.0)
    self.assertLess(recipe.optimized_l_over_r, math.sqrt(math.pi))

  def test_pin_material_map_uses_central_mask_unless_pin_is_moderator(self):
    recipe = optimized_pin_recipe_7()
    uo2 = material_map_for_pin(0, recipe=recipe)
    moderator = material_map_for_pin(6, recipe=recipe)
    self.assertEqual(len(uo2), 49)
    self.assertEqual(uo2.count(0), 13)
    self.assertEqual(uo2.count(6), 36)
    self.assertEqual(set(moderator), {6})

  def test_assembly_and_core_expansions_have_expected_sizes(self):
    recipe = optimized_pin_recipe_7()
    assembly_materials = expanded_materials(
      ASSEMBLY_PIN_MAPS[0],
      pins_per_side=ASSEMBLY_DIMENSION,
      recipe=recipe,
    )
    self.assertEqual(
      len(assembly_materials),
      (ASSEMBLY_DIMENSION * recipe.number) ** 2,
    )
    pin_map = core_pin_map(
      CORE_ASSEMBLY_MAP,
      ASSEMBLY_PIN_MAPS,
      assemblies_per_side=CORE_DIMENSION,
      pins_per_assembly_side=ASSEMBLY_DIMENSION,
    )
    self.assertEqual(
      len(pin_map),
      (CORE_DIMENSION * ASSEMBLY_DIMENSION) ** 2,
    )
    core_edges = repeated_edges(
      recipe.widths,
      CORE_DIMENSION * ASSEMBLY_DIMENSION,
    )
    self.assertEqual(len(core_edges), CORE_DIMENSION * ASSEMBLY_DIMENSION * recipe.number + 1)
    self.assertAlmostEqual(core_edges[-1], CORE_DIMENSION * ASSEMBLY_DIMENSION * PIN_PITCH)

  def test_center_sampled_pin_recipe_conserves_area_with_density_corrections(self):
    recipe = center_sampled_mesh_recipe([0], pins_per_side=1, number=3)
    corrections = {
      entry.material: entry
      for entry in recipe.material_area_corrections
    }

    self.assertTrue(recipe.conservative)
    self.assertEqual(center_sampled_material_map_for_pin(0, number=3).count(0), 5)
    for entry in corrections.values():
      corrected_area = entry.assigned_area * entry.density_correction
      self.assertAlmostEqual(corrected_area, entry.target_area)
    self.assertAlmostEqual(corrections[0].target_area, math.pi * PIN_RADIUS * PIN_RADIUS)
    self.assertAlmostEqual(
      corrections[6].target_area,
      PIN_PITCH * PIN_PITCH - math.pi * PIN_RADIUS * PIN_RADIUS,
    )

  def test_center_sampled_recipe_flags_missing_materials(self):
    recipe = center_sampled_mesh_recipe([0], pins_per_side=1, number=1)
    corrections = {
      entry.material: entry
      for entry in recipe.material_area_corrections
    }

    self.assertFalse(recipe.conservative)
    self.assertIsNone(corrections[6].density_correction)
    self.assertGreater(corrections[6].target_area, 0.0)
    self.assertEqual(corrections[6].assigned_area, 0.0)

  def test_density_corrected_materials_scale_cross_sections_only(self):
    base = get_materials()
    corrected = get_density_corrected_materials({0: 1.25, 6: 0.75})

    self.assertAlmostEqual(corrected.sigma_t(0, 0), 1.25 * base.sigma_t(0, 0))
    self.assertAlmostEqual(corrected.sigma_f(0, 6), 1.25 * base.sigma_f(0, 6))
    self.assertAlmostEqual(corrected.sigma_s(6, 5, 6), 0.75 * base.sigma_s(6, 5, 6))
    self.assertAlmostEqual(corrected.nu(0, 0), base.nu(0, 0))
    self.assertAlmostEqual(corrected.chi(0, 0), base.chi(0, 0))


if __name__ == "__main__":
  unittest.main()
