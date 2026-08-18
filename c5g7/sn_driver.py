"""Artifact-writing SN/DD drivers for C5G7 pin, assembly, and core cases."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import detran as dt

try:
  from .assemblies_c5g7 import get_assemblies
  from .core_c5g7 import get_core
  from .csg_c5g7 import (
    ASSEMBLY_PIN_MAPS,
    CORE_ASSEMBLY_MAP,
    PIN_NAMES,
  )
  from .material_c5g7 import get_density_corrected_materials, get_materials
  from .pins_c5g7 import get_pins
  from .reference import assembly_eigenvalues, eigenvalue, pin_eigenvalues
  from .sn_mesh import (
    ASSEMBLY_DIMENSION,
    CORE_DIMENSION,
    build_center_sampled_problem_mesh,
    build_lattice_mesh,
    build_pin_mesh,
    core_pin_map,
    mesh_audit,
    optimized_pin_recipe_7,
  )
except ImportError:
  from assemblies_c5g7 import get_assemblies
  from core_c5g7 import get_core
  from csg_c5g7 import ASSEMBLY_PIN_MAPS, CORE_ASSEMBLY_MAP, PIN_NAMES
  from material_c5g7 import get_density_corrected_materials, get_materials
  from pins_c5g7 import get_pins
  from reference import assembly_eigenvalues, eigenvalue, pin_eigenvalues
  from sn_mesh import (
    ASSEMBLY_DIMENSION,
    CORE_DIMENSION,
    build_center_sampled_problem_mesh,
    build_lattice_mesh,
    build_pin_mesh,
    core_pin_map,
    mesh_audit,
    optimized_pin_recipe_7,
  )


SN_OPTION_NAMES = (
  "index",
  "mesh_variant",
  "pin_mesh",
  "density_correction",
  "equation",
  "quadrature",
  "polar",
  "azimuth",
  "eigen_backend",
  "acceleration",
  "linear_backend",
  "linear_tolerance",
  "linear_iters",
  "eigen_tolerance",
  "eigen_iters",
  "inner_iters",
  "outer_method",
  "outer_iters",
  "outer_tolerance",
  "cmfd_eigen_algorithm",
  "cmfd_coarse_mesh_level",
  "cmfd_relaxation",
  "monitor_level",
  "print_interval",
)

DEFAULT_SN_OPTIONS = {
  "index": 0,
  "mesh_variant": "optimized",
  "pin_mesh": 7,
  "density_correction": True,
  "equation": "dd",
  "quadrature": "asqr-asdr",
  "polar": 3,
  "azimuth": 6,
  "eigen_backend": "slepc-gd",
  "acceleration": "mgcmdsa",
  "linear_backend": "petsc",
  "linear_tolerance": 1.0e-9,
  "linear_iters": 5000,
  "eigen_tolerance": 1.0e-6,
  "eigen_iters": 1000,
  "inner_iters": 1,
  "outer_method": "gmres",
  "outer_iters": 1000,
  "outer_tolerance": 1.0e-8,
  "cmfd_eigen_algorithm": "power",
  "cmfd_coarse_mesh_level": None,
  "cmfd_relaxation": 0.7,
  "monitor_level": 1,
  "print_interval": 2,
}

CASE_SN_OPTIONS = {
  "pin": {
    "polar": 5,
    "azimuth": 10,
    "quadrature": "u-u",
    "outer_iters": 100,
    "eigen_tolerance": 1.0e-8,
  },
  "assembly": {
    "polar": 5,
    "azimuth": 10,
    "quadrature": "u-u",
    "outer_iters": 400,
    "eigen_tolerance": 1.0e-8,
  },
  "core": {},
}

C5G7_CORE_BOUNDARY = ("reflect", "vacuum", "reflect", "vacuum")


def reference_for_case(case: str, index: int):
  if case == "pin":
    return pin_eigenvalues().get(index)
  if case == "assembly":
    return assembly_eigenvalues().get(index)
  if case == "core":
    return ("C5G7 2-D core reference", eigenvalue())
  raise ValueError("case must be pin, assembly, or core")


def boundary_for_case(case: str):
  if case in ("pin", "assembly"):
    return "reflect"
  if case == "core":
    return C5G7_CORE_BOUNDARY
  raise ValueError("case must be pin, assembly, or core")


def _legacy_mesh(case: str, index: int, pin_mesh: int):
  if case == "pin":
    pins = get_pins(pin_mesh, True)
    if index < 0 or index >= len(pins):
      raise ValueError("pin index must be between 0 and 6")
    return pins[index].mesh()
  if case == "assembly":
    assemblies = get_assemblies(pin_mesh, True)
    if index < 0 or index >= len(assemblies):
      raise ValueError("assembly index must be between 0 and 2")
    return assemblies[index].mesh()
  if case == "core":
    if index not in (None, 0):
      raise ValueError("core index must be 0")
    return get_core(pin_mesh, True).mesh()
  raise ValueError("case must be pin, assembly, or core")


def _optimized_mesh(case: str, index: int):
  recipe = optimized_pin_recipe_7()
  if case == "pin":
    if index < 0 or index > 6:
      raise ValueError("pin index must be between 0 and 6")
    return build_pin_mesh(dt, index, recipe), recipe
  if case == "assembly":
    if index < 0 or index >= len(ASSEMBLY_PIN_MAPS):
      raise ValueError("assembly index must be between 0 and 2")
    return (
      build_lattice_mesh(
        dt,
        ASSEMBLY_PIN_MAPS[index],
        pins_per_side=ASSEMBLY_DIMENSION,
        recipe=recipe,
      ),
      recipe,
    )
  if case == "core":
    if index not in (None, 0):
      raise ValueError("core index must be 0")
    pin_map = core_pin_map(
      CORE_ASSEMBLY_MAP,
      ASSEMBLY_PIN_MAPS,
      assemblies_per_side=CORE_DIMENSION,
      pins_per_assembly_side=ASSEMBLY_DIMENSION,
    )
    return (
      build_lattice_mesh(
        dt,
        pin_map,
        pins_per_side=CORE_DIMENSION * ASSEMBLY_DIMENSION,
        recipe=recipe,
      ),
      recipe,
    )
  raise ValueError("case must be pin, assembly, or core")


def _centered_mesh(case: str, index: int, pin_mesh: int):
  if case == "pin":
    if index < 0 or index > 6:
      raise ValueError("pin index must be between 0 and 6")
    return build_center_sampled_problem_mesh(
      dt,
      [index],
      pins_per_side=1,
      number=pin_mesh,
    )
  if case == "assembly":
    if index < 0 or index >= len(ASSEMBLY_PIN_MAPS):
      raise ValueError("assembly index must be between 0 and 2")
    return build_center_sampled_problem_mesh(
      dt,
      ASSEMBLY_PIN_MAPS[index],
      pins_per_side=ASSEMBLY_DIMENSION,
      number=pin_mesh,
    )
  if case == "core":
    if index not in (None, 0):
      raise ValueError("core index must be 0")
    pin_map = core_pin_map(
      CORE_ASSEMBLY_MAP,
      ASSEMBLY_PIN_MAPS,
      assemblies_per_side=CORE_DIMENSION,
      pins_per_assembly_side=ASSEMBLY_DIMENSION,
    )
    return build_center_sampled_problem_mesh(
      dt,
      pin_map,
      pins_per_side=CORE_DIMENSION * ASSEMBLY_DIMENSION,
      number=pin_mesh,
    )
  raise ValueError("case must be pin, assembly, or core")


def get_sn_case(case: str, index: int, *, mesh_variant: str, pin_mesh: int):
  """Return mesh, mesh metadata, and reference data for an SN case."""

  if mesh_variant == "legacy":
    mesh = _legacy_mesh(case, index, pin_mesh)
    recipe = {
      "name": "legacy_pin_cell_meshify",
      "pin_mesh": pin_mesh,
      "volume_preserving": bool(pin_mesh >= 7),
      "source": "PinCell.meshify(number, True)",
    }
  elif mesh_variant == "optimized":
    if pin_mesh != 7:
      raise ValueError("the automated optimized mesh is currently implemented for 7x7")
    mesh, recipe_obj = _optimized_mesh(case, index)
    recipe = recipe_obj.to_dict()
  elif mesh_variant == "centered":
    mesh, recipe_obj = _centered_mesh(case, index, pin_mesh)
    recipe = recipe_obj.to_dict()
  else:
    raise ValueError("mesh_variant must be legacy, optimized, or centered")
  return mesh, recipe, reference_for_case(case, index)


def material_for_mesh_recipe(recipe: dict, *, apply_density_correction=True):
  corrections = recipe.get("density_correction_by_material") or {}
  if not corrections or not apply_density_correction:
    return get_materials(), {
      "name": "C5G7",
      "density_corrected": False,
      "density_correction_applied": False,
      "density_correction_available": bool(corrections),
      "conservative": True,
    }
  finite = {
    int(material): float(factor)
    for material, factor in corrections.items()
    if factor is not None
  }
  return get_density_corrected_materials(finite), {
    "name": "C5G7 density-corrected",
    "density_corrected": True,
    "density_correction_applied": True,
    "density_correction_available": True,
    "correction_scope": recipe.get("correction_scope"),
    "conservative": bool(recipe.get("conservative_density_correction")),
    "density_correction_by_material": corrections,
  }


def _linear_solver(options: dict):
  linear = options["linear_backend"]
  common = {
    "rtol": options["linear_tolerance"],
    "atol": options["linear_tolerance"],
    "max_iters": options["linear_iters"],
    "monitor_level": 0,
    "pc_type": "petsc_pc",
    "petsc_pc_type": "ilu",
  }
  if linear == "petsc":
    return dt.LinearSolver.petsc(**common)
  if linear == "gmres":
    common.pop("petsc_pc_type")
    common["pc_type"] = "ilu0"
    return dt.LinearSolver.gmres(restart=30, **common)
  raise ValueError("linear_backend must be petsc or gmres")


def _solver_database(options: dict):
  if options["eigen_backend"] == "cmfd":
    eigen_algorithm = options["cmfd_eigen_algorithm"]
  elif options["eigen_backend"].startswith("slepc"):
    eigen_algorithm = "slepc"
  else:
    eigen_algorithm = "gd"
  return dt.SolverDatabase(
    linear=_linear_solver(options),
    eigen_algorithm=eigen_algorithm,
    eigen_tolerance=options["eigen_tolerance"],
    eigen_max_iters=options["eigen_iters"],
    eigen_monitor_level=options["monitor_level"],
    options={
      "eigen_solver_slepc_type": (
        "gd" if options["eigen_backend"] == "slepc-gd" else "krylovschur"
      )
    } if options["eigen_backend"].startswith("slepc") else {},
  )


def _eigen_solver(options: dict):
  backend = options["eigen_backend"]
  common = {
    "tolerance": options["eigen_tolerance"],
    "max_iters": options["eigen_iters"],
    "monitor_level": options["monitor_level"],
    "print_interval": max(1, options["print_interval"]),
    "inner": _linear_solver(options),
  }
  if backend == "slepc-gd":
    return dt.EigenSolver.slepc(method="arnoldi", slepc_type="gd", **common)
  if backend == "slepc-krylovschur":
    return dt.EigenSolver.slepc(
      method="arnoldi",
      slepc_type="krylovschur",
      **common,
    )
  if backend == "gd":
    return dt.EigenSolver(method="GD", algorithm="gd", **common)
  if backend == "cmfd":
    return dt.EigenSolver(
      method="cmfd",
      algorithm=options["cmfd_eigen_algorithm"],
      **common,
    )
  if backend == "power":
    return dt.EigenSolver.power_iteration(**common)
  raise ValueError("unsupported eigen backend")


def _outer_solver(options: dict):
  method = options["outer_method"]
  if method == "source":
    return dt.IterationSolver.source_iteration(
      tolerance=options["outer_tolerance"],
      max_iters=options["outer_iters"],
      print_level=0,
      print_interval=max(1, options["print_interval"]),
    )
  if method == "gauss-seidel":
    return dt.IterationSolver.gauss_seidel(
      tolerance=options["outer_tolerance"],
      max_iters=options["outer_iters"],
      print_level=0,
      print_interval=max(1, options["print_interval"]),
    )
  if method != "gmres":
    raise ValueError("outer_method must be gmres, gauss-seidel, or source")
  acceleration = options["acceleration"]
  if acceleration == "cmfd":
    return dt.IterationSolver(
      "CMFD",
      tolerance=options["outer_tolerance"],
      max_iters=options["outer_iters"],
      print_level=0,
      print_interval=max(1, options["print_interval"]),
    )
  pc_type = "mgcmdsa" if acceleration == "mgcmdsa" else "none"
  return dt.IterationSolver.gmres(
    tolerance=options["outer_tolerance"],
    max_iters=options["outer_iters"],
    print_level=0,
    print_interval=max(1, options["print_interval"]),
    krylov_group_cutoff=0,
    pc_type=pc_type,
  )


def build_sn_problem(mesh, material, *, case: str, options: dict):
  if options["eigen_backend"] == "cmfd" and options["acceleration"] != "cmfd":
    raise ValueError("CMFD eigen backend requires --acceleration cmfd")
  callow_db = _solver_database(options)
  overrides = {
    "inner_print_level": 0,
    "outer_print_level": 0,
  }
  if options["acceleration"] == "mgcmdsa":
    overrides.update({
      "outer_pc_type": "mgcmdsa",
      "eigen_solver_pc_type": "mgcmdsa",
      "mgpc_coarse_mesh_level": options["pin_mesh"],
      "mgpc_condensation_option": 1,
      "mgpc_cmdsa_use_smoothing": 1,
      "mgpc_cmdsa_smoothing_iters": 2,
      "mgpc_cmdsa_smoothing_relax": 0.8,
    })
  if options["acceleration"] == "cmfd":
    level = options["cmfd_coarse_mesh_level"]
    if level is None:
      level = options["pin_mesh"]
    overrides.update({
      "mgpc_coarse_mesh_level": level,
      "cmfd_relaxation": options["cmfd_relaxation"],
    })

  return dt.EigenvalueProblem(
    mesh,
    material,
    equation=options["equation"],
    boundary=boundary_for_case(case),
    inner=dt.IterationSolver.gmres(
      max_iters=options["inner_iters"],
      tolerance=options["outer_tolerance"],
      print_level=0,
    ),
    outer=_outer_solver(options),
    solver=_eigen_solver(options),
    quadrature=dt.Quadrature.product(
      polar=options["polar"],
      azimuth=options["azimuth"],
      type=options["quadrature"],
    ),
    solver_databases={
      "inner_solver_db": callow_db,
      "inner_pc_db": callow_db,
      "outer_solver_db": callow_db,
      "outer_pc_db": callow_db,
      "eigen_solver_db": callow_db,
      "eigen_solver_pc_db": callow_db,
    },
    input_overrides=overrides,
  )


def _mesh_summary(mesh) -> dict:
  return {
    "cells": int(mesh.number_cells()),
    "cells_x": int(mesh.number_cells_x()),
    "cells_y": int(mesh.number_cells_y()),
    "width_x": float(mesh.total_width_x()),
    "width_y": float(mesh.total_width_y()),
  }


def _case_name(case: str, index: int) -> str:
  if case == "pin":
    return PIN_NAMES.get(index, f"pin {index}")
  if case == "assembly":
    return ("UO2 assembly", "MOX assembly", "reflector assembly")[index]
  return "C5G7 quarter core"


def summarize_sn_case(case: str, index: int, *, options: dict, audit_only: bool = False):
  mesh, recipe, reference = get_sn_case(
    case,
    index,
    mesh_variant=options["mesh_variant"],
    pin_mesh=options["pin_mesh"],
  )
  material, material_model = material_for_mesh_recipe(
    recipe,
    apply_density_correction=options["density_correction"],
  )
  result = {
    "case": case,
    "index": index,
    "name": _case_name(case, index),
    "method": "native_sn",
    "geometry_model": (
      "cartesian_center_sampled_density_corrected_pin_mesh"
      if recipe.get("center_sampled") and material_model["density_correction_applied"]
      else "cartesian_center_sampled_uncorrected_pin_mesh"
      if recipe.get("center_sampled")
      else "cartesian_volume_preserving_pin_mesh"
    ),
    "mesh": _mesh_summary(mesh),
    "pin_mesh_recipe": recipe,
    "material_model": material_model,
    "discretization": {
      "equation": options["equation"],
      "quadrature": options["quadrature"],
      "polar_per_octant": options["polar"],
      "azimuth_per_octant": options["azimuth"],
      "boundary": boundary_for_case(case),
    },
    "solver": {
      "eigen_backend": options["eigen_backend"],
      "acceleration": options["acceleration"],
      "linear_backend": options["linear_backend"],
      "linear_tolerance": options["linear_tolerance"],
      "linear_iters": options["linear_iters"],
      "eigen_tolerance": options["eigen_tolerance"],
      "eigen_iters": options["eigen_iters"],
      "inner_iters": options["inner_iters"],
      "outer_method": options["outer_method"],
      "outer_iters": options["outer_iters"],
      "outer_tolerance": options["outer_tolerance"],
      "cmfd_eigen_algorithm": options["cmfd_eigen_algorithm"],
      "cmfd_coarse_mesh_level": options["cmfd_coarse_mesh_level"],
      "cmfd_relaxation": options["cmfd_relaxation"],
    },
    "reference": {},
    "results": {},
  }
  if reference is not None:
    result["reference"] = {
      "name": reference[0],
      "keff": reference[1],
      "source": "C5G7 benchmark reference",
    }
  if audit_only:
    result["results"]["status"] = "audit_only"
    return result

  problem = build_sn_problem(mesh, material, case=case, options=options)
  start = time.time()
  dt.Manager.initialize([sys.argv[0]])
  solve_result = problem.solve()
  keff = float(solve_result.keff)
  sweeps = None
  if hasattr(solve_result.solver, "number_sweeps"):
    sweeps = int(solve_result.solver.number_sweeps())
  elapsed = time.time() - start
  ref_keff = result["reference"].get("keff")
  result["results"] = {
    "status": "solved",
    "keff": keff,
    "finite": math.isfinite(keff),
    "elapsed_s": elapsed,
    "manager_finalized": False,
  }
  if ref_keff is not None:
    result["results"]["delta_keff"] = keff - ref_keff
    result["results"]["delta_pcm"] = 1.0e5 * (keff - ref_keff)
  if sweeps is not None:
    result["results"]["sweeps"] = sweeps
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


def latex_result_row(result: dict) -> str:
  results = result.get("results", {})
  reference = result.get("reference", {})
  values = [
    result["name"],
    result["pin_mesh_recipe"]["name"],
    result["mesh"]["cells"],
    result["discretization"]["equation"],
    result["discretization"]["quadrature"],
    result["discretization"]["polar_per_octant"],
    result["discretization"]["azimuth_per_octant"],
    result["solver"]["eigen_backend"],
    result["solver"]["acceleration"],
    _format_float(reference.get("keff"), 9),
    _format_float(results.get("keff"), 9),
    _format_float(results.get("delta_pcm"), 4),
    _format_float(results.get("elapsed_s"), 4),
  ]
  return " & ".join(_latex_escape(value) for value in values)


def write_case_artifacts(result: dict, output_dir: Path, stem: str) -> dict:
  output_dir.mkdir(parents=True, exist_ok=True)
  json_path = output_dir / f"{stem}.json"
  row_path = output_dir / f"{stem}_row.tex"
  json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
  row_path.write_text(latex_result_row(result) + "\n")
  return {"json": str(json_path), "latex_row": str(row_path)}


def add_sn_arguments(
  parser: argparse.ArgumentParser,
  *,
  include_case: bool = False,
  default_case: str = "core",
  default_index: int = 0,
) -> argparse.ArgumentParser:
  if include_case:
    parser.add_argument("--case", choices=("pin", "assembly", "core"), default=default_case)
  parser.add_argument("--index", type=int, default=default_index)
  parser.add_argument(
    "--mesh-variant",
    choices=("optimized", "legacy", "centered"),
    default="optimized",
  )
  parser.add_argument("--pin-mesh", type=int, default=7)
  parser.add_argument("--density-correction", dest="density_correction",
                      action="store_true", default=None)
  parser.add_argument("--no-density-correction", dest="density_correction",
                      action="store_false")
  parser.add_argument("--equation", choices=("dd", "sc"), default="dd")
  parser.add_argument("--quadrature", default=None)
  parser.add_argument("--polar", type=int, default=None)
  parser.add_argument("--azimuth", type=int, default=None)
  parser.add_argument(
    "--eigen-backend",
    choices=("slepc-gd", "slepc-krylovschur", "gd", "cmfd", "power"),
    default=None,
  )
  parser.add_argument("--acceleration", choices=("mgcmdsa", "cmfd", "none"), default=None)
  parser.add_argument("--linear-backend", choices=("petsc", "gmres"), default=None)
  parser.add_argument("--linear-tolerance", type=float, default=None)
  parser.add_argument("--linear-iters", type=int, default=None)
  parser.add_argument("--eigen-tolerance", type=float, default=None)
  parser.add_argument("--eigen-iters", type=int, default=None)
  parser.add_argument("--inner-iters", type=int, default=None)
  parser.add_argument(
    "--outer-method",
    choices=("gmres", "gauss-seidel", "source"),
    default=None,
  )
  parser.add_argument("--outer-iters", type=int, default=None)
  parser.add_argument("--outer-tolerance", type=float, default=None)
  parser.add_argument("--cmfd-eigen-algorithm", choices=("power", "gd"), default=None)
  parser.add_argument("--cmfd-coarse-mesh-level", type=int, default=None)
  parser.add_argument("--cmfd-relaxation", type=float, default=None)
  parser.add_argument("--monitor-level", type=int, default=None)
  parser.add_argument("--print-interval", type=int, default=None)
  parser.add_argument("--audit-only", action="store_true")
  parser.add_argument("--output-dir", default=None)
  parser.add_argument("--stem", default=None)
  return parser


def options_from_args(case: str, args) -> dict:
  options = dict(DEFAULT_SN_OPTIONS)
  options.update(CASE_SN_OPTIONS.get(case, {}))
  for name in SN_OPTION_NAMES:
    value = getattr(args, name, None)
    if value is not None:
      options[name] = value
  return options


def run_from_args(args, *, case: str | None = None) -> dict:
  selected_case = case or args.case
  options = options_from_args(selected_case, args)
  result = summarize_sn_case(
    selected_case,
    options["index"],
    options=options,
    audit_only=args.audit_only,
  )
  stem = args.stem or f"sn_{selected_case}_{options['index']}_{options['mesh_variant']}"
  output_dir = (
    Path(args.output_dir)
    if args.output_dir
    else Path(__file__).resolve().parent / "results" / "sn"
  )
  paths = write_case_artifacts(result, output_dir, stem)
  result["artifacts"] = paths
  return result


def print_summary(result: dict) -> None:
  print(json.dumps(result, indent=2, sort_keys=True))


def main(default_case: str = "core", default_index: int = 0, default_stem: str | None = None) -> None:
  parser = add_sn_arguments(
    argparse.ArgumentParser(description="Run or audit a C5G7 SN/DD case."),
    include_case=True,
    default_case=default_case,
    default_index=default_index,
  )
  args = parser.parse_args()
  if default_stem is not None and args.stem is None:
    args.stem = default_stem
  print_summary(run_from_args(args, case=args.case))


if __name__ == "__main__":
  main()
