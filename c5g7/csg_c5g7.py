"""True 2-D CSG geometry scaffold for the C5G7 benchmark.

The existing C5G7 drivers use a Cartesian pin mesh.  This module builds the
same pin, assembly, and quarter-core layouts as analytic CSG regions.  Pin
subdivision is controlled by radial divisions in the fuel/central material,
radial divisions in the moderator, and azimuthal sectors.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass
import inspect
import json
import math
from pathlib import Path
import time
from typing import Iterable

import numpy as np

import detran as dt

try:
  from .assemblies_c5g7 import get_assemblies
  from .core_c5g7 import get_core
  from .material_c5g7 import get_materials
  from .pins_c5g7 import get_pins
  from .reference import assembly_eigenvalues, eigenvalue, pin_eigenvalues, pin_powers
  from .sn_mesh import (
    build_center_sampled_problem_mesh,
    build_lattice_mesh,
    build_pin_mesh,
    core_pin_map,
    optimized_pin_recipe_7,
  )
except ImportError:
  from assemblies_c5g7 import get_assemblies
  from core_c5g7 import get_core
  from material_c5g7 import get_materials
  from pins_c5g7 import get_pins
  from reference import assembly_eigenvalues, eigenvalue, pin_eigenvalues, pin_powers
  from sn_mesh import (
    build_center_sampled_problem_mesh,
    build_lattice_mesh,
    build_pin_mesh,
    core_pin_map,
    optimized_pin_recipe_7,
  )


PIN_PITCH = 1.26
PIN_RADIUS = 0.54
ASSEMBLY_DIMENSION = 17
ASSEMBLY_PITCH = ASSEMBLY_DIMENSION * PIN_PITCH
CORE_DIMENSION = 3
CORE_PITCH = CORE_DIMENSION * ASSEMBLY_PITCH

PIN_NAMES = {
  0: "UO2 fuel",
  1: "MOX 4.3% fuel",
  2: "MOX 7.0% fuel",
  3: "MOX 8.7% fuel",
  4: "guide tube",
  5: "fission chamber",
  6: "moderator",
}

ASSEMBLY_NAMES = {
  0: "UO2 assembly",
  1: "MOX assembly",
  2: "reflector assembly",
}

CORE_ASSEMBLY_MAP = (
  0, 1, 2,
  1, 0, 2,
  2, 2, 2,
)

C5G7_CORE_BOUNDARY = ("reflect", "vacuum", "reflect", "vacuum")

UO2_PIN_MAP = (
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,4,0,0,4,0,0,4,0,0,0,0,0,
  0,0,0,4,0,0,0,0,0,0,0,0,0,4,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,4,0,0,4,0,0,4,0,0,4,0,0,4,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,4,0,0,4,0,0,5,0,0,4,0,0,4,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,4,0,0,4,0,0,4,0,0,4,0,0,4,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,4,0,0,0,0,0,0,0,0,0,4,0,0,0,
  0,0,0,0,0,4,0,0,4,0,0,4,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
)

MOX_PIN_MAP = (
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
  1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,
  1,2,2,2,2,4,2,2,4,2,2,4,2,2,2,2,1,
  1,2,2,4,2,3,3,3,3,3,3,3,2,4,2,2,1,
  1,2,2,2,3,3,3,3,3,3,3,3,3,2,2,2,1,
  1,2,4,3,3,4,3,3,4,3,3,4,3,3,4,2,1,
  1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
  1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
  1,2,4,3,3,4,3,3,5,3,3,4,3,3,4,2,1,
  1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
  1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
  1,2,4,3,3,4,3,3,4,3,3,4,3,3,4,2,1,
  1,2,2,2,3,3,3,3,3,3,3,3,3,2,2,2,1,
  1,2,2,4,2,3,3,3,3,3,3,3,2,4,2,2,1,
  1,2,2,2,2,4,2,2,4,2,2,4,2,2,2,2,1,
  1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
)

REFLECTOR_PIN_MAP = (6,) * (ASSEMBLY_DIMENSION * ASSEMBLY_DIMENSION)
ASSEMBLY_PIN_MAPS = (UO2_PIN_MAP, MOX_PIN_MAP, REFLECTOR_PIN_MAP)


@dataclass(frozen=True)
class PinDiscretization:
  """CSG pin subdivision controls."""

  fuel_radial_divisions: int = 1
  moderator_radial_divisions: int = 1
  azimuthal_slices: int = 8
  azimuthal_offset_degrees: float = 0.0

  def validate(self) -> None:
    if self.fuel_radial_divisions < 1:
      raise ValueError("fuel_radial_divisions must be positive")
    if self.moderator_radial_divisions < 1:
      raise ValueError("moderator_radial_divisions must be positive")
    if self.azimuthal_slices < 3:
      raise ValueError("azimuthal_slices must be at least 3")
    if 360.0 / self.azimuthal_slices >= 180.0:
      raise ValueError("azimuthal sector width must be less than 180 degrees")


@dataclass(frozen=True)
class RegionRecord:
  """Metadata matching the region insertion order used by CSGModel2D."""

  region: int
  material: int
  material_name: str
  pin: int
  pin_name: str
  pin_i: int
  pin_j: int
  assembly: int | None = None
  assembly_name: str | None = None
  assembly_i: int | None = None
  assembly_j: int | None = None
  zone: str = "fuel"
  radial_index: int = 0
  sector_index: int = 0
  theta0_degrees: float = 0.0
  theta1_degrees: float = 0.0


@dataclass(frozen=True)
class CSGC5G7Case:
  """A built C5G7 CSG model plus Python-side region metadata."""

  name: str
  case: str
  index: int
  model: object
  discretization: PinDiscretization
  regions: tuple[RegionRecord, ...]


def _radial_edges(start: float, stop: float, divisions: int) -> list[float]:
  return [
    start + (stop - start) * i / divisions
    for i in range(divisions + 1)
  ]


def _sector_bounds(discretization: PinDiscretization) -> list[tuple[float, float]]:
  width = 360.0 / discretization.azimuthal_slices
  start = discretization.azimuthal_offset_degrees
  return [
    (math.radians(start + i * width), math.radians(start + (i + 1) * width))
    for i in range(discretization.azimuthal_slices)
  ]


def _add_region(records: list[RegionRecord], **kwargs) -> None:
  records.append(RegionRecord(region=len(records), **kwargs))


def _pin_material_is_homogeneous(pin: int) -> bool:
  return pin == 6


def _add_csg_pin(
  model,
  records: list[RegionRecord],
  pin: int,
  x0: float,
  y0: float,
  discretization: PinDiscretization,
  *,
  pin_i: int,
  pin_j: int,
  assembly: int | None = None,
  assembly_i: int | None = None,
  assembly_j: int | None = None,
) -> None:
  cx = x0 + 0.5 * PIN_PITCH
  cy = y0 + 0.5 * PIN_PITCH
  corner_radius = math.sqrt(2.0) * 0.5 * PIN_PITCH
  sector_bounds = _sector_bounds(discretization)
  assembly_name = ASSEMBLY_NAMES.get(assembly) if assembly is not None else None

  if not _pin_material_is_homogeneous(pin):
    fuel_edges = _radial_edges(0.0, PIN_RADIUS, discretization.fuel_radial_divisions)
    for radial_index, (r0, r1) in enumerate(zip(fuel_edges[:-1], fuel_edges[1:])):
      for sector_index, (theta0, theta1) in enumerate(sector_bounds):
        model.add_radial_sector(pin, cx, cy, r0, r1, theta0, theta1)
        _add_region(
          records,
          material=pin,
          material_name=PIN_NAMES[pin],
          pin=pin,
          pin_name=PIN_NAMES[pin],
          pin_i=pin_i,
          pin_j=pin_j,
          assembly=assembly,
          assembly_name=assembly_name,
          assembly_i=assembly_i,
          assembly_j=assembly_j,
          zone="fuel",
          radial_index=radial_index,
          sector_index=sector_index,
          theta0_degrees=math.degrees(theta0),
          theta1_degrees=math.degrees(theta1),
        )
    moderator_start = PIN_RADIUS
  else:
    moderator_start = 0.0

  moderator_edges = _radial_edges(
    moderator_start,
    corner_radius,
    discretization.moderator_radial_divisions,
  )
  for radial_index, (r0, r1) in enumerate(zip(moderator_edges[:-1], moderator_edges[1:])):
    for sector_index, (theta0, theta1) in enumerate(sector_bounds):
      model.add_rectangle_radial_sector(
        6,
        x0,
        x0 + PIN_PITCH,
        y0,
        y0 + PIN_PITCH,
        cx,
        cy,
        r0,
        r1,
        theta0,
        theta1,
      )
      _add_region(
        records,
        material=6,
        material_name=PIN_NAMES[6],
        pin=pin,
        pin_name=PIN_NAMES[pin],
        pin_i=pin_i,
        pin_j=pin_j,
        assembly=assembly,
        assembly_name=assembly_name,
        assembly_i=assembly_i,
        assembly_j=assembly_j,
        zone="moderator",
        radial_index=radial_index,
        sector_index=sector_index,
        theta0_degrees=math.degrees(theta0),
        theta1_degrees=math.degrees(theta1),
      )


def build_pin(index: int = 0, discretization: PinDiscretization | None = None) -> CSGC5G7Case:
  """Build one reflected C5G7 pin cell as true 2-D CSG."""

  discretization = discretization or PinDiscretization()
  discretization.validate()
  if index < 0 or index > 6:
    raise ValueError("pin index must be between 0 and 6")

  model = dt.CSGModel2D(PIN_PITCH, PIN_PITCH)
  records: list[RegionRecord] = []
  _add_csg_pin(
    model,
    records,
    index,
    0.0,
    0.0,
    discretization,
    pin_i=0,
    pin_j=0,
  )
  return CSGC5G7Case(
    name=PIN_NAMES[index],
    case="pin",
    index=index,
    model=model,
    discretization=discretization,
    regions=tuple(records),
  )


def build_assembly(
  index: int = 0,
  discretization: PinDiscretization | None = None,
) -> CSGC5G7Case:
  """Build one C5G7 assembly as true 2-D CSG."""

  discretization = discretization or PinDiscretization()
  discretization.validate()
  if index < 0 or index >= len(ASSEMBLY_PIN_MAPS):
    raise ValueError("assembly index must be between 0 and 2")

  model = dt.CSGModel2D(ASSEMBLY_PITCH, ASSEMBLY_PITCH)
  records: list[RegionRecord] = []
  pin_map = ASSEMBLY_PIN_MAPS[index]
  for j in range(ASSEMBLY_DIMENSION):
    for i in range(ASSEMBLY_DIMENSION):
      pin = pin_map[j * ASSEMBLY_DIMENSION + i]
      _add_csg_pin(
        model,
        records,
        pin,
        i * PIN_PITCH,
        j * PIN_PITCH,
        discretization,
        pin_i=i,
        pin_j=j,
        assembly=index,
      )
  return CSGC5G7Case(
    name=ASSEMBLY_NAMES[index],
    case="assembly",
    index=index,
    model=model,
    discretization=discretization,
    regions=tuple(records),
  )


def build_core(discretization: PinDiscretization | None = None) -> CSGC5G7Case:
  """Build the 3x3 C5G7 quarter-core as true 2-D CSG."""

  discretization = discretization or PinDiscretization()
  discretization.validate()
  model = dt.CSGModel2D(CORE_PITCH, CORE_PITCH)
  records: list[RegionRecord] = []
  for aj in range(CORE_DIMENSION):
    for ai in range(CORE_DIMENSION):
      assembly = CORE_ASSEMBLY_MAP[aj * CORE_DIMENSION + ai]
      pin_map = ASSEMBLY_PIN_MAPS[assembly]
      assembly_x0 = ai * ASSEMBLY_PITCH
      assembly_y0 = aj * ASSEMBLY_PITCH
      for pj in range(ASSEMBLY_DIMENSION):
        for pi in range(ASSEMBLY_DIMENSION):
          pin = pin_map[pj * ASSEMBLY_DIMENSION + pi]
          _add_csg_pin(
            model,
            records,
            pin,
            assembly_x0 + pi * PIN_PITCH,
            assembly_y0 + pj * PIN_PITCH,
            discretization,
            pin_i=pi,
            pin_j=pj,
            assembly=assembly,
            assembly_i=ai,
            assembly_j=aj,
          )
  return CSGC5G7Case(
    name="C5G7 quarter core",
    case="core",
    index=0,
    model=model,
    discretization=discretization,
    regions=tuple(records),
  )


def build_case(
  case: str,
  index: int = 0,
  discretization: PinDiscretization | None = None,
) -> CSGC5G7Case:
  if case == "pin":
    return build_pin(index, discretization)
  if case == "assembly":
    return build_assembly(index, discretization)
  if case == "core":
    if index not in (None, 0):
      raise ValueError("core index must be 0")
    return build_core(discretization)
  raise ValueError("case must be pin, assembly, or core")


def reference_for_case(case: str, index: int) -> dict:
  if case == "pin":
    ref = pin_eigenvalues().get(index)
  elif case == "assembly":
    ref = assembly_eigenvalues().get(index)
  elif case == "core":
    ref = ("C5G7 2-D core reference", eigenvalue())
  else:
    ref = None
  if ref is None:
    return {}
  return {"name": ref[0], "keff": ref[1], "source": "C5G7 benchmark reference"}


def boundary_for_case(case: str):
  """Return the 2-D boundary condition for a C5G7 CSG transport case."""

  if case in ("pin", "assembly"):
    return "reflect"
  if case == "core":
    return C5G7_CORE_BOUNDARY
  raise ValueError("case must be pin, assembly, or core")


def _material_counts(records: Iterable[RegionRecord]) -> dict[str, int]:
  counts: dict[str, int] = {}
  for record in records:
    key = str(record.material)
    counts[key] = counts.get(key, 0) + 1
  return counts


def estimate_region_areas(product, number_regions: int) -> list[float]:
  """Estimate CSG region areas by length times track spatial weight."""

  direction_weights = product.directions.weight
  track_direction = product.tracks.direction_index
  track_weight = product.tracks.spatial_weight
  segment_track = product.segments.track_index
  weighted_area = track_weight[segment_track] * direction_weights[track_direction[segment_track]]
  total_direction_weight = float(np.sum(direction_weights))
  areas = np.zeros(number_regions, dtype=float)
  for region in np.unique(product.segments.region_id):
    mask = product.segments.region_id == region
    areas[int(region)] = float(
      np.sum(product.segments.length[mask] * weighted_area[mask])
      / total_direction_weight
    )
  return [float(value) for value in areas]


def material_areas(records: Iterable[RegionRecord], region_areas: list[float]) -> dict[str, float]:
  areas: dict[str, float] = {}
  for record in records:
    key = str(record.material)
    areas[key] = areas.get(key, 0.0) + float(region_areas[record.region])
  return areas


def _pin_map_for_case(case: str, index: int) -> tuple[list[int], int]:
  """Return the C5G7 pin map and number of pins per side for a case."""

  if case == "pin":
    if index < 0 or index > 6:
      raise ValueError("pin index must be between 0 and 6")
    return [index], 1
  if case == "assembly":
    if index < 0 or index >= len(ASSEMBLY_PIN_MAPS):
      raise ValueError("assembly index must be between 0 and 2")
    return list(ASSEMBLY_PIN_MAPS[index]), ASSEMBLY_DIMENSION
  if case == "core":
    if index not in (None, 0):
      raise ValueError("core index must be 0")
    return core_pin_map(
      CORE_ASSEMBLY_MAP,
      ASSEMBLY_PIN_MAPS,
      assemblies_per_side=CORE_DIMENSION,
      pins_per_assembly_side=ASSEMBLY_DIMENSION,
    ), CORE_DIMENSION * ASSEMBLY_DIMENSION
  raise ValueError("case must be pin, assembly, or core")


def _legacy_overlay_mesh(case: str, index: int, pin_mesh: int):
  if case == "pin":
    return get_pins(pin_mesh, True)[index].mesh()
  if case == "assembly":
    return get_assemblies(pin_mesh, True)[index].mesh()
  if case == "core":
    return get_core(pin_mesh, True).mesh()
  raise ValueError("case must be pin, assembly, or core")


def build_csg_overlay_mesh(case: str, index: int, *, variant: str, pin_mesh: int):
  """Return a Cartesian overlay mesh and metadata for CSG FSR refinement."""

  variant = str(variant)
  pin_mesh = int(pin_mesh)
  if variant == "none":
    return None, {
      "enabled": False,
      "variant": "none",
      "note": "No Cartesian overlay; each base CSG region is one FSR.",
    }
  if variant == "legacy":
    mesh = _legacy_overlay_mesh(case, index, pin_mesh)
    return mesh, {
      "enabled": True,
      "variant": "legacy",
      "pin_mesh": pin_mesh,
      "source": "PinCell.meshify(number, True)",
      "cells": int(mesh.number_cells()),
      "width_x": float(mesh.total_width_x()),
      "width_y": float(mesh.total_width_y()),
      "material_map_note": (
        "Only the Cartesian mesh planes are used for CSG overlay splitting; "
        "materials come from the analytic CSG region crossed by each track."
      ),
    }

  pin_map, pins_per_side = _pin_map_for_case(case, index)
  if variant == "optimized":
    if pin_mesh != 7:
      raise ValueError("optimized CSG overlay mesh is currently implemented for 7x7")
    recipe = optimized_pin_recipe_7()
    if case == "pin":
      mesh = build_pin_mesh(dt, index, recipe)
    else:
      mesh = build_lattice_mesh(
        dt,
        pin_map,
        pins_per_side=pins_per_side,
        recipe=recipe,
      )
    metadata = recipe.to_dict()
    metadata.update({
      "enabled": True,
      "variant": "optimized",
      "pin_mesh": pin_mesh,
      "cells": int(mesh.number_cells()),
      "width_x": float(mesh.total_width_x()),
      "width_y": float(mesh.total_width_y()),
      "material_map_note": (
        "Only the optimized Cartesian mesh planes are used for CSG overlay "
        "splitting; materials remain true analytic CSG materials."
      ),
    })
    return mesh, metadata

  if variant == "centered":
    mesh, recipe = build_center_sampled_problem_mesh(
      dt,
      pin_map,
      pins_per_side=pins_per_side,
      number=pin_mesh,
    )
    metadata = recipe.to_dict()
    metadata.update({
      "enabled": True,
      "variant": "centered",
      "pin_mesh": pin_mesh,
      "cells": int(mesh.number_cells()),
      "width_x": float(mesh.total_width_x()),
      "width_y": float(mesh.total_width_y()),
      "material_map_note": (
        "Only the centered Cartesian mesh planes are used for CSG overlay "
        "splitting; materials remain true analytic CSG materials, and density "
        "corrections or center-sampled materials are not used in the true-CSG solve."
      ),
    })
    return mesh, metadata

  raise ValueError("csg overlay mesh variant must be none, legacy, optimized, or centered")


def run_mesh_moc_reference(case: str, index: int, args) -> dict:
  """Run the existing mesh-MOC backend for comparison."""

  try:
    from .moc_driver import run_moc_case
  except ImportError:
    from moc_driver import run_moc_case

  summary = run_moc_case(
    case,
    index=index,
    pin_mesh=args.mesh_pin_mesh,
    polar=args.moc_polar,
    azimuth=args.moc_azimuth,
    quadrature=args.moc_quadrature,
    spacing=args.moc_spacing,
    equation=args.moc_equation,
    inner_tolerance=args.moc_inner_tolerance,
    outer_tolerance=args.moc_outer_tolerance,
    eigen_tolerance=args.moc_eigen_tolerance,
    inner_iters=args.moc_inner_iters,
    outer_iters=args.moc_outer_iters,
    eigen_iters=args.moc_eigen_iters,
    monitor_level=args.moc_monitor_level,
    print_interval=args.moc_print_interval,
    plot_dir=None,
  )
  return {
    "backend": "mesh_moc",
    "note": "Existing staircase-mesh MOC backend for comparison.",
    **summary,
  }


def _track_collection_counts(tracks) -> dict:
  """Return counts for a raw CSG TrackCollection exposed through pydetran."""

  tables = dt.csg_track_collection_tables(tracks)
  return {
    "directions": int(len(tables.direction_key0)),
    "tracks": int(len(tables.track_direction_index)),
    "segments": int(len(tables.segment_track_index)),
  }


def _cell_fission_rate(material, volumes, fluxes, cell: int, record: RegionRecord) -> float:
  """Return sigma_f-weighted fission rate for one synthetic FSR cell."""

  volume = float(volumes[cell])
  rate = 0.0
  for group, flux in enumerate(fluxes):
    rate += (
      float(material.sigma_f(record.material, group))
      * float(flux[cell])
      * volume
    )
  return rate


def _normalized_pin_power_errors(
  calculated: np.ndarray,
  reference: np.ndarray,
) -> dict:
  """Normalize calculated pin powers to the reference active-pin sum."""

  active = reference > 0.0
  if calculated.shape != reference.shape:
    raise ValueError("calculated and reference pin-power grids must match")
  active_calculated_sum = float(np.sum(calculated[active]))
  if active_calculated_sum <= 0.0:
    raise ValueError("calculated active-pin power sum must be positive")

  normalized = calculated * (float(np.sum(reference[active])) / active_calculated_sum)
  relative_error = np.full(reference.shape, np.nan, dtype=float)
  relative_error[active] = (
    (normalized[active] - reference[active]) / reference[active] * 100.0
  )
  abs_error = np.abs(relative_error[active])
  active_indices = np.argwhere(active)
  max_local = int(np.argmax(abs_error))
  max_j, max_i = (int(value) for value in active_indices[max_local])
  return {
    "normalization": "calculated active-pin fission rates scaled to reference active-pin sum",
    "active_pin_count": int(np.count_nonzero(active)),
    "reference_sum": float(np.sum(reference[active])),
    "calculated_raw_sum": active_calculated_sum,
    "rms_percent": float(np.sqrt(np.mean(relative_error[active] ** 2))),
    "mean_abs_percent": float(np.mean(abs_error)),
    "max_abs_percent": float(abs_error[max_local]),
    "max_abs_pin": {
      "i": max_i,
      "j": max_j,
      "reference": float(reference[max_j, max_i]),
      "calculated": float(normalized[max_j, max_i]),
      "relative_error_percent": float(relative_error[max_j, max_i]),
    },
    "calculated_normalized": normalized.tolist(),
    "reference": reference.tolist(),
    "relative_error_percent": relative_error.tolist(),
  }


def core_pin_power_summary(csg_case: CSGC5G7Case, material, region_mesh, solve_result) -> dict | None:
  """Compute normalized quarter-core pin powers and errors for a CSG solve."""

  if csg_case.case != "core":
    return None

  pins_per_side = CORE_DIMENSION * ASSEMBLY_DIMENSION
  rates = np.zeros((pins_per_side, pins_per_side), dtype=float)
  records_by_region = {record.region: record for record in csg_case.regions}
  volumes = np.asarray(
    [float(region_mesh.volume(region)) for region in range(region_mesh.number_cells())],
    dtype=float,
  )
  fluxes = [
    np.asarray(solve_result.flux(group), dtype=float)
    for group in range(int(material.number_groups()))
  ]
  if region_mesh.mesh_map_exists("CSG_REGION"):
    csg_regions = [int(value) for value in region_mesh.mesh_map("CSG_REGION")]
  else:
    csg_regions = list(range(int(region_mesh.number_cells())))

  for cell, csg_region in enumerate(csg_regions):
    record = records_by_region[int(csg_region)]
    if record.assembly_i is None or record.assembly_j is None:
      continue
    i = record.assembly_i * ASSEMBLY_DIMENSION + record.pin_i
    j = record.assembly_j * ASSEMBLY_DIMENSION + record.pin_j
    rates[j, i] += _cell_fission_rate(material, volumes, fluxes, cell, record)

  reference = np.asarray(pin_powers(), dtype=float).reshape(
    pins_per_side,
    pins_per_side,
  )
  summary = _normalized_pin_power_errors(rates, reference)
  summary.update({
    "case": "core",
    "rate_model": (
      "sum_g sigma_f(material,g) * phi_g(FSR cell) * FSR volume, "
      "accumulated through CSG region metadata"
    ),
    "grid_shape": [pins_per_side, pins_per_side],
    "calculated_raw_fission_rate": rates.tolist(),
  })
  return summary


def run_native_csg_moc(csg_case: CSGC5G7Case, args) -> dict:
  """Run native CSG-region MOC through the Detran eigenvalue manager."""

  reference = reference_for_case(csg_case.case, csg_case.index)
  boundary = boundary_for_case(csg_case.case)
  material = get_materials()
  overlay_mesh, overlay_metadata = build_csg_overlay_mesh(
    csg_case.case,
    csg_case.index,
    variant=args.csg_overlay_mesh_variant,
    pin_mesh=args.csg_overlay_pin_mesh,
  )
  problem_kwargs = {}
  if overlay_mesh is not None:
    parameters = inspect.signature(dt.CSGEigenvalueProblem.__init__).parameters
    if "overlay_mesh" not in parameters:
      raise RuntimeError(
        "this pydetran build does not expose CSGEigenvalueProblem overlay_mesh; "
        "use the newer build-pydetran wrapper or run without a CSG overlay"
      )
    problem_kwargs["overlay_mesh"] = overlay_mesh
  problem = dt.CSGEigenvalueProblem(
    csg_case.model,
    material,
    equation=args.moc_equation,
    boundary=boundary,
    inner=dt.IterationSolver.source_iteration(
      tolerance=args.moc_inner_tolerance,
      max_iters=args.moc_inner_iters,
      print_level=0,
    ),
    outer=dt.IterationSolver.gauss_seidel(
      tolerance=args.moc_outer_tolerance,
      max_iters=args.moc_outer_iters,
      print_level=0,
    ),
    solver=dt.EigenSolver.power_iteration(
      tolerance=args.moc_eigen_tolerance,
      max_iters=args.moc_eigen_iters,
      monitor_level=args.moc_monitor_level,
      print_interval=max(1, args.moc_print_interval),
    ),
    quadrature=dt.Quadrature.product(
      polar=args.moc_polar,
      azimuth=args.moc_azimuth,
      type=args.moc_quadrature,
    ),
    input_overrides={
      "tracker_maximum_spacing": args.moc_spacing,
      "tracker_minimum_number": 1,
      "tracker_spatial_quad_type": "uniform",
      "inner_print_level": 0,
      "outer_print_level": 0,
    },
    **problem_kwargs,
  )

  start = time.time()
  solve_result = problem.solve()
  elapsed = time.time() - start
  keff = float(solve_result.keff)
  summary = {
    "backend": "native_csg_moc",
    "status": "completed",
    "fsr_model": (
      "true_csg_regions_split_by_cartesian_overlay"
      if overlay_mesh is not None else
      "one_fsr_per_true_csg_region"
    ),
    "volume_model": (
      "True CSG tracks are split at Cartesian overlay planes. Each FSR is a "
      "non-empty (CSG region, overlay cell) intersection; FSR volumes are "
      "estimated from the CSG overlay track product and scaled to the analytic "
      "domain area."
      if overlay_mesh is not None else
      "One synthetic mesh cell per CSG region; region volumes are estimated "
      "from the CSG track product and scaled to the analytic domain area."
    ),
    "base_csg_regions": int(csg_case.model.number_regions),
    "overlay_mesh": overlay_metadata,
    "boundary": boundary,
    "cells": int(problem.region_mesh.number_cells()),
    "polar": args.moc_polar,
    "azimuth": args.moc_azimuth,
    "quadrature": args.moc_quadrature,
    "spacing": args.moc_spacing,
    "equation": args.moc_equation,
    "keff": keff,
    "finite": math.isfinite(keff),
    "elapsed": elapsed,
  }
  if hasattr(solve_result.solver, "number_sweeps"):
    summary["sweeps"] = int(solve_result.solver.number_sweeps())
  if problem.tracks is not None:
    summary.update(_track_collection_counts(problem.tracks))
  overlay_product = getattr(problem, "overlay_product", None)
  if overlay_product is not None:
    summary["overlay_product"] = {
      "fsrs": int(len(overlay_product.csg_region)),
      "unique_csg_regions": int(len(set(int(value) for value in overlay_product.csg_region))),
      "unique_overlay_cells": int(len(set(int(value) for value in overlay_product.overlay_cell))),
      "has_region_mesh_maps": all(
        problem.region_mesh.mesh_map_exists(key)
        for key in ("CSG_REGION", "OVERLAY_CELL", "OVERLAY_I", "OVERLAY_J")
      ),
    }
  if reference:
    summary["reference_name"] = reference["name"]
    summary["reference_keff"] = reference["keff"]
    summary["delta_keff"] = keff - reference["keff"]
  pin_power = core_pin_power_summary(
    csg_case,
    material,
    problem.region_mesh,
    solve_result,
  )
  if pin_power is not None:
    summary["pin_power"] = pin_power
  return summary


def summarize_case(csg_case: CSGC5G7Case, *, args=None) -> dict:
  """Return a JSON-serializable benchmark summary for one CSG case."""

  model = csg_case.model
  result = {
    "case": csg_case.case,
    "index": csg_case.index,
    "name": csg_case.name,
    "created_by": "c5g7.csg_c5g7",
    "geometry_backend": "native_2d_csg",
    "transport_backend": None,
    "transport_status": (
      "Native CSG geometry/tracking is built here. Add --solve-csg-moc to "
      "run the region-track MOC eigenvalue bridge."
    ),
    "discretization": asdict(csg_case.discretization),
    "width_x": float(model.width_x),
    "width_y": float(model.width_y),
    "number_regions": int(model.number_regions),
    "region_metadata_count": len(csg_case.regions),
    "material_region_counts": _material_counts(csg_case.regions),
    "reference": reference_for_case(csg_case.case, csg_case.index),
  }

  if args is not None and args.track:
    angles = [float(value) for value in args.track_angles.split(",") if value.strip()]
    t0 = time.perf_counter()
    product = model.track(
      angles=angles,
      spacing=args.track_spacing,
      minimum_tracks=args.track_minimum,
    )
    elapsed = time.perf_counter() - t0
    result["track_audit"] = {
      "angles_degrees": angles,
      "spacing": args.track_spacing,
      "minimum_tracks": args.track_minimum,
      "tracks": int(len(product.tracks)),
      "segments": int(len(product.segments)),
      "elapsed_s": elapsed,
    }

  if args is not None and args.solve_csg_moc:
    result["transport"] = run_native_csg_moc(csg_case, args)
    result["transport_backend"] = result["transport"]["backend"]
    result["transport_status"] = "Native CSG-region MOC solve completed."

  if args is not None and args.solve_mesh_moc:
    comparison = run_mesh_moc_reference(csg_case.case, csg_case.index, args)
    if result.get("transport_backend"):
      result["mesh_moc_comparison"] = comparison
    else:
      result["transport"] = comparison
      result["transport_backend"] = comparison["backend"]
      result["transport_status"] = "Mesh-MOC comparison solve completed."

  return result


def _latex_escape(value) -> str:
  text = "" if value is None else str(value)
  return (
    text.replace("\\", r"\textbackslash{}")
    .replace("&", r"\&")
    .replace("%", r"\%")
    .replace("_", r"\_")
    .replace("#", r"\#")
  )


def _format_float(value, digits: int = 6) -> str:
  if value is None:
    return "--"
  return f"{float(value):.{digits}g}"


def _backend_label(value) -> str:
  labels = {
    "native_csg_moc": "CSG-MOC",
    "mesh_moc": "mesh-MOC",
  }
  return labels.get(value, value or "--")


def latex_result_row(result: dict) -> str:
  reference = result.get("reference", {})
  transport = result.get("transport") or {}
  keff = transport.get("keff")
  reference_keff = transport.get("reference_keff", reference.get("keff"))
  delta_pcm = None
  if keff is not None and reference_keff is not None:
    delta_pcm = 1.0e5 * (float(keff) - float(reference_keff))
  track = result.get("track_audit", {})
  disc = result["discretization"]
  values = [
    result["name"],
    disc["fuel_radial_divisions"],
    disc["moderator_radial_divisions"],
    disc["azimuthal_slices"],
    result["number_regions"],
    track.get("tracks", "--"),
    track.get("segments", "--"),
    _backend_label(transport.get("backend", result.get("transport_backend"))),
    _format_float(reference_keff, 9),
    _format_float(keff, 9),
    _format_float(delta_pcm, 4),
  ]
  return " & ".join(_latex_escape(value) for value in values)


def write_case_artifacts(result: dict, output_dir: Path, stem: str) -> dict:
  output_dir.mkdir(parents=True, exist_ok=True)
  json_path = output_dir / f"{stem}.json"
  row_path = output_dir / f"{stem}_row.tex"
  json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
  row_path.write_text(latex_result_row(result) + "\n")
  return {"json": str(json_path), "latex_row": str(row_path)}


def _console_result(result: dict) -> dict:
  """Return a compact result view for stdout."""

  view = copy.deepcopy(result)
  pin_power = (view.get("transport") or {}).get("pin_power")
  if pin_power:
    for key in (
      "calculated_normalized",
      "reference",
      "relative_error_percent",
      "calculated_raw_fission_rate",
    ):
      pin_power.pop(key, None)
  return view


def add_csg_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
  parser.add_argument("--case", choices=("pin", "assembly", "core"))
  parser.add_argument("--index", type=int, default=0)
  parser.add_argument("--fuel-radial-divisions", type=int, default=1)
  parser.add_argument("--moderator-radial-divisions", type=int, default=1)
  parser.add_argument("--azimuthal-slices", type=int, default=8)
  parser.add_argument("--azimuthal-offset-degrees", type=float, default=0.0)
  parser.add_argument("--output-dir", default=None)
  parser.add_argument("--stem", default=None)
  parser.add_argument("--track", action="store_true")
  parser.add_argument("--track-angles", default="0,45,90,135")
  parser.add_argument("--track-spacing", type=float, default=2.0)
  parser.add_argument("--track-minimum", type=int, default=1)
  parser.add_argument("--solve-csg-moc", action="store_true")
  parser.add_argument("--csg-overlay-mesh-variant",
                      choices=("none", "legacy", "optimized", "centered"),
                      default="none",
                      help="Cartesian mesh used only to split true CSG tracks into overlay FSRs")
  parser.add_argument("--csg-overlay-pin-mesh", type=int, default=7,
                      help="pin-cell subdivisions for centered/legacy CSG overlay meshes")
  parser.add_argument("--solve-mesh-moc", action="store_true")
  parser.add_argument("--mesh-pin-mesh", type=int, default=7)
  parser.add_argument("--moc-polar", type=int, default=1)
  parser.add_argument("--moc-azimuth", type=int, default=4)
  parser.add_argument("--moc-quadrature", default="u-u")
  parser.add_argument("--moc-spacing", type=float, default=0.25)
  parser.add_argument("--moc-equation", choices=("scmoc", "ddmoc"), default="scmoc")
  parser.add_argument("--moc-inner-tolerance", type=float, default=1.0e-5)
  parser.add_argument("--moc-outer-tolerance", type=float, default=1.0e-5)
  parser.add_argument("--moc-eigen-tolerance", type=float, default=1.0e-4)
  parser.add_argument("--moc-inner-iters", type=int, default=80)
  parser.add_argument("--moc-outer-iters", type=int, default=80)
  parser.add_argument("--moc-eigen-iters", type=int, default=60)
  parser.add_argument("--moc-monitor-level", type=int, default=1)
  parser.add_argument("--moc-print-interval", type=int, default=10)
  return parser


def main(default_case: str = "core", default_index: int = 0, default_stem: str | None = None) -> None:
  parser = add_csg_arguments(
    argparse.ArgumentParser(description="Build and audit a C5G7 2-D CSG case.")
  )
  parser.set_defaults(case=default_case, index=default_index)
  args = parser.parse_args()
  discretization = PinDiscretization(
    fuel_radial_divisions=args.fuel_radial_divisions,
    moderator_radial_divisions=args.moderator_radial_divisions,
    azimuthal_slices=args.azimuthal_slices,
    azimuthal_offset_degrees=args.azimuthal_offset_degrees,
  )
  csg_case = build_case(args.case, args.index, discretization)
  result = summarize_case(csg_case, args=args)
  stem = args.stem or default_stem or f"csg_{args.case}_{args.index}"
  output_dir = Path(args.output_dir) if args.output_dir else Path(__file__).resolve().parent / "csg_results"
  paths = write_case_artifacts(result, output_dir, stem)
  result["artifacts"] = paths
  print(json.dumps(_console_result(result), indent=2, sort_keys=True))


if __name__ == "__main__":
  main()
