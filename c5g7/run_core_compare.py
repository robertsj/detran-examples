# pyexamples/c5g7/run_core_compare.py
#
# Compare native MOC and SN on the same 2-D C5G7 quarter-core mesh.

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import detran as dt

try:
  from .core_c5g7 import get_core
  from .material_c5g7 import get_materials
  from .moc_driver import (
    C5G7_CORE_BOUNDARY,
    build_moc_problem,
    parse_plot_groups,
  )
  from .reference import eigenvalue
except ImportError:
  from core_c5g7 import get_core
  from material_c5g7 import get_materials
  from moc_driver import (
    C5G7_CORE_BOUNDARY,
    build_moc_problem,
    parse_plot_groups,
  )
  from reference import eigenvalue


DEFAULT_OPTIONS = {
  "pin_mesh": 1,
  "polar": 1,
  "azimuth": 1,
  "quadrature": "u-u",
  "spacing": 2.0,
  "inner_tolerance": 1.0e-3,
  "outer_tolerance": 1.0e-3,
  "eigen_tolerance": 1.0e-3,
  "inner_iters": 20,
  "outer_iters": 20,
  "eigen_iters": 8,
  "monitor_level": 0,
  "print_interval": 1,
  "moc_equation": "scmoc",
  "sn_equation": "dd",
  "edit_map": "PINS",
  "groups": None,
  "normalization": "global-volume",
  "npz_output": None,
}


def build_sn_problem(
  mesh,
  material,
  *,
  polar=1,
  azimuth=1,
  quadrature="u-u",
  inner_tolerance=1.0e-3,
  outer_tolerance=1.0e-3,
  eigen_tolerance=1.0e-3,
  inner_iters=20,
  outer_iters=20,
  eigen_iters=8,
  monitor_level=0,
  print_interval=1,
  boundary=C5G7_CORE_BOUNDARY,
  equation="dd",
):
  """Return an SN eigenproblem matched to the MOC mesh and scalar unknowns."""
  return dt.EigenvalueProblem(
    mesh,
    material,
    equation=equation,
    boundary=boundary,
    inner=dt.IterationSolver.source_iteration(
      tolerance=inner_tolerance,
      max_iters=inner_iters,
      print_level=0,
    ),
    outer=dt.IterationSolver.gauss_seidel(
      tolerance=outer_tolerance,
      max_iters=outer_iters,
      print_level=0,
    ),
    solver=dt.EigenSolver.power_iteration(
      tolerance=eigen_tolerance,
      max_iters=eigen_iters,
      monitor_level=monitor_level,
      print_interval=max(1, print_interval),
    ),
    quadrature=dt.Quadrature.product(
      polar=polar,
      azimuth=azimuth,
      type=quadrature,
    ),
    input_overrides={
      "inner_print_level": 0,
      "outer_print_level": 0,
    },
  )


def _mesh_volumes(mesh):
  return np.fromiter(
    (mesh.volume(i) for i in range(int(mesh.number_cells()))),
    dtype=float,
    count=int(mesh.number_cells()),
  )


def _flux_matrix(result, groups):
  return np.vstack([
    np.asarray(result.flux(group), dtype=float)
    for group in groups
  ])


def _normalization_scale(flux, volumes, mode):
  if mode == "none":
    return 1.0
  total = float(np.sum(flux * volumes[None, :]))
  if not math.isfinite(total) or abs(total) < 1.0e-300:
    return 1.0
  return total


def _normalized_flux(flux, volumes, mode):
  return flux / _normalization_scale(flux, volumes, mode)


def _labels(mesh, edit_map):
  if edit_map is None or str(edit_map).lower() == "cell":
    return "cell", np.arange(int(mesh.number_cells()), dtype=np.int64)
  key = str(edit_map)
  if not mesh.mesh_map_exists(key):
    raise ValueError("mesh map '{0}' does not exist".format(key))
  return key, np.asarray(mesh.mesh_map(key), dtype=np.int64)


def _volume_average_by_label(flux, labels, volumes):
  unique, inverse = np.unique(labels, return_inverse=True)
  label_volume = np.bincount(inverse, weights=volumes, minlength=unique.size)
  averages = []
  for group_flux in flux:
    weighted = np.bincount(
      inverse,
      weights=group_flux * volumes,
      minlength=unique.size,
    )
    averages.append(np.divide(
      weighted,
      label_volume,
      out=np.zeros_like(weighted),
      where=label_volume > 0.0,
    ))
  return unique, label_volume, np.vstack(averages)


def _group_stats(moc_flux, sn_flux, volumes, groups):
  stats = []
  volume = float(np.sum(volumes))
  for row, group in enumerate(groups):
    diff = moc_flux[row] - sn_flux[row]
    abs_diff = np.abs(diff)
    weighted_l2 = math.sqrt(float(np.sum(diff * diff * volumes) / volume))
    sn_l2 = math.sqrt(float(np.sum(sn_flux[row] * sn_flux[row] * volumes) / volume))
    stats.append({
      "group": int(group),
      "mean_abs": float(np.sum(abs_diff * volumes) / volume),
      "max_abs": float(np.max(abs_diff)),
      "l2": weighted_l2,
      "relative_l2": float(weighted_l2 / sn_l2) if sn_l2 > 0.0 else None,
    })
  return stats


def _solver_summary(name, result, elapsed):
  summary = {
    "name": name,
    "keff": float(result.keff),
    "finite": bool(math.isfinite(result.keff)),
    "elapsed": float(elapsed),
  }
  if hasattr(result.solver, "number_sweeps"):
    summary["sweeps"] = int(result.solver.number_sweeps())
  return summary


def run(**options):
  opts = dict(DEFAULT_OPTIONS)
  opts.update(options)

  material = get_materials()
  core = get_core(int(opts["pin_mesh"]), True)
  mesh = core.mesh()
  groups = opts["groups"]
  if groups is None:
    groups = tuple(range(int(material.number_groups())))
  else:
    groups = tuple(int(group) for group in groups)

  boundary = C5G7_CORE_BOUNDARY
  common = {
    "polar": opts["polar"],
    "azimuth": opts["azimuth"],
    "quadrature": opts["quadrature"],
    "inner_tolerance": opts["inner_tolerance"],
    "outer_tolerance": opts["outer_tolerance"],
    "eigen_tolerance": opts["eigen_tolerance"],
    "inner_iters": opts["inner_iters"],
    "outer_iters": opts["outer_iters"],
    "eigen_iters": opts["eigen_iters"],
    "monitor_level": opts["monitor_level"],
    "print_interval": opts["print_interval"],
    "boundary": boundary,
  }
  moc_problem = build_moc_problem(
    mesh,
    material,
    spacing=opts["spacing"],
    equation=opts["moc_equation"],
    **common,
  )
  sn_problem = build_sn_problem(
    mesh,
    material,
    equation=opts["sn_equation"],
    **common,
  )

  start = time.time()
  moc_result = moc_problem.solve()
  moc_elapsed = time.time() - start

  start = time.time()
  sn_result = sn_problem.solve()
  sn_elapsed = time.time() - start

  volumes = _mesh_volumes(mesh)
  moc_flux = _flux_matrix(moc_result, groups)
  sn_flux = _flux_matrix(sn_result, groups)
  if moc_flux.shape != sn_flux.shape:
    raise RuntimeError(
      "MOC/SN flux shapes differ: {0} vs {1}".format(
        moc_flux.shape,
        sn_flux.shape,
      )
    )

  normalization = opts["normalization"]
  moc_cmp = _normalized_flux(moc_flux, volumes, normalization)
  sn_cmp = _normalized_flux(sn_flux, volumes, normalization)

  cell_label_name, cell_labels = _labels(mesh, "cell")
  edit_label_name, edit_labels = _labels(mesh, opts["edit_map"])
  edit_ids, edit_volumes, moc_edit = _volume_average_by_label(
    moc_cmp,
    edit_labels,
    volumes,
  )
  _, _, sn_edit = _volume_average_by_label(sn_cmp, edit_labels, volumes)

  summary = {
    "case": "C5G7 core",
    "reference_keff": float(eigenvalue()),
    "pin_mesh": int(opts["pin_mesh"]),
    "cells": int(mesh.number_cells()),
    "number_cells_x": int(mesh.number_cells_x()),
    "number_cells_y": int(mesh.number_cells_y()),
    "groups": [int(group) for group in groups],
    "boundary": tuple(boundary),
    "polar": int(opts["polar"]),
    "azimuth": int(opts["azimuth"]),
    "quadrature": opts["quadrature"],
    "spacing": float(opts["spacing"]),
    "moc_equation": opts["moc_equation"],
    "sn_equation": opts["sn_equation"],
    "normalization": normalization,
    "cell_scalar_unknowns": int(mesh.number_cells()) * len(groups),
    "edit_map": edit_label_name,
    "edit_regions": int(edit_ids.size),
    "edit_scalar_unknowns": int(edit_ids.size) * len(groups),
    "moc": _solver_summary("MOC", moc_result, moc_elapsed),
    "sn": _solver_summary("SN", sn_result, sn_elapsed),
    "delta_keff_moc_minus_sn": float(moc_result.keff - sn_result.keff),
    "delta_keff_moc_minus_reference": float(moc_result.keff - eigenvalue()),
    "delta_keff_sn_minus_reference": float(sn_result.keff - eigenvalue()),
    "cell_comparison": _group_stats(moc_cmp, sn_cmp, volumes, groups),
    "edit_comparison": _group_stats(moc_edit, sn_edit, edit_volumes, groups),
  }

  if opts["npz_output"]:
    path = Path(opts["npz_output"])
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
      path,
      groups=np.asarray(groups, dtype=np.int64),
      volumes=volumes,
      cell_labels=cell_labels,
      edit_labels=edit_labels,
      edit_ids=edit_ids,
      moc_flux=moc_flux,
      sn_flux=sn_flux,
      moc_flux_normalized=moc_cmp,
      sn_flux_normalized=sn_cmp,
      moc_edit_flux=moc_edit,
      sn_edit_flux=sn_edit,
    )
    summary["npz_output"] = str(path)

  return summary


def _format_float(value):
  if value is None:
    return "n/a"
  return "{:.8e}".format(float(value))


def print_summary(summary):
  print("case: {0}".format(summary["case"]))
  print("boundary(W,E,S,N): {0}".format(summary["boundary"]))
  print("mesh: {0} x {1} = {2} cells".format(
    summary["number_cells_x"],
    summary["number_cells_y"],
    summary["cells"],
  ))
  print("groups: {0}".format(summary["groups"]))
  print("cell scalar unknowns: {0}".format(summary["cell_scalar_unknowns"]))
  print("edit map: {0} ({1} regions, {2} scalar unknowns)".format(
    summary["edit_map"],
    summary["edit_regions"],
    summary["edit_scalar_unknowns"],
  ))
  print("quadrature: {0}, polar={1}, azimuth={2}".format(
    summary["quadrature"],
    summary["polar"],
    summary["azimuth"],
  ))
  print("MOC equation: {0}".format(summary["moc_equation"]))
  print("SN equation: {0}".format(summary["sn_equation"]))
  print("normalization: {0}".format(summary["normalization"]))
  print("reference keff: {:.9f}".format(summary["reference_keff"]))
  print("moc keff: {:.9f} elapsed={:.3f}s sweeps={}".format(
    summary["moc"]["keff"],
    summary["moc"]["elapsed"],
    summary["moc"].get("sweeps", "n/a"),
  ))
  print("sn  keff: {:.9f} elapsed={:.3f}s sweeps={}".format(
    summary["sn"]["keff"],
    summary["sn"]["elapsed"],
    summary["sn"].get("sweeps", "n/a"),
  ))
  print("delta keff (MOC-SN): {:.8e}".format(summary["delta_keff_moc_minus_sn"]))
  print("")
  _print_stats("cell normalized flux comparison (MOC - SN)",
               summary["cell_comparison"])
  if summary["edit_map"] != "cell":
    print("")
    _print_stats(
      "{0} volume-averaged normalized flux comparison (MOC - SN)".format(
        summary["edit_map"],
      ),
      summary["edit_comparison"],
    )
  if "npz_output" in summary:
    print("")
    print("npz_output: {0}".format(summary["npz_output"]))


def _print_stats(title, stats):
  print(title)
  print("group  rel_l2       max_abs      mean_abs")
  for item in stats:
    print("{group:>5d}  {rel}  {max_abs}  {mean_abs}".format(
      group=item["group"],
      rel=_format_float(item["relative_l2"]),
      max_abs=_format_float(item["max_abs"]),
      mean_abs=_format_float(item["mean_abs"]),
    ))


def parse_args():
  parser = argparse.ArgumentParser(
    description="Compare native 2-D MOC and SN on the full C5G7 core mesh.",
  )
  parser.add_argument("--pin-mesh", type=int, default=DEFAULT_OPTIONS["pin_mesh"],
                      help="pin-cell mesh subdivisions per dimension")
  parser.add_argument("--polar", type=int, default=DEFAULT_OPTIONS["polar"])
  parser.add_argument("--azimuth", type=int, default=DEFAULT_OPTIONS["azimuth"])
  parser.add_argument("--quadrature", default=DEFAULT_OPTIONS["quadrature"])
  parser.add_argument("--spacing", type=float, default=DEFAULT_OPTIONS["spacing"],
                      help="MOC maximum average track spacing")
  parser.add_argument("--inner-tolerance", type=float,
                      default=DEFAULT_OPTIONS["inner_tolerance"])
  parser.add_argument("--outer-tolerance", type=float,
                      default=DEFAULT_OPTIONS["outer_tolerance"])
  parser.add_argument("--eigen-tolerance", type=float,
                      default=DEFAULT_OPTIONS["eigen_tolerance"])
  parser.add_argument("--inner-iters", type=int,
                      default=DEFAULT_OPTIONS["inner_iters"])
  parser.add_argument("--outer-iters", type=int,
                      default=DEFAULT_OPTIONS["outer_iters"])
  parser.add_argument("--eigen-iters", type=int,
                      default=DEFAULT_OPTIONS["eigen_iters"])
  parser.add_argument("--monitor-level", type=int,
                      default=DEFAULT_OPTIONS["monitor_level"])
  parser.add_argument("--print-interval", type=int,
                      default=DEFAULT_OPTIONS["print_interval"])
  parser.add_argument("--moc-equation", choices=("scmoc", "ddmoc"),
                      default=DEFAULT_OPTIONS["moc_equation"])
  parser.add_argument("--sn-equation", choices=("dd", "sc", "sd"),
                      default=DEFAULT_OPTIONS["sn_equation"])
  parser.add_argument("--edit-map", default=DEFAULT_OPTIONS["edit_map"],
                      help="mesh map to reduce over, or 'cell'")
  parser.add_argument("--groups", type=parse_plot_groups, default=None,
                      help="zero-based group list/ranges, e.g. 0,2-4 or all")
  parser.add_argument("--normalization", choices=("global-volume", "none"),
                      default=DEFAULT_OPTIONS["normalization"])
  parser.add_argument("--npz-output", default=None,
                      help="optional NPZ output containing raw flux arrays")
  parser.add_argument("--json", action="store_true",
                      help="print the full summary as JSON")
  return parser.parse_args()


def _groups_from_args(value):
  if value == "all":
    return None
  return value


if __name__ == "__main__":
  args = parse_args()
  summary = run(
    pin_mesh=args.pin_mesh,
    polar=args.polar,
    azimuth=args.azimuth,
    quadrature=args.quadrature,
    spacing=args.spacing,
    inner_tolerance=args.inner_tolerance,
    outer_tolerance=args.outer_tolerance,
    eigen_tolerance=args.eigen_tolerance,
    inner_iters=args.inner_iters,
    outer_iters=args.outer_iters,
    eigen_iters=args.eigen_iters,
    monitor_level=args.monitor_level,
    print_interval=args.print_interval,
    moc_equation=args.moc_equation,
    sn_equation=args.sn_equation,
    edit_map=args.edit_map,
    groups=_groups_from_args(args.groups),
    normalization=args.normalization,
    npz_output=args.npz_output,
  )
  if args.json:
    print(json.dumps(summary, indent=2, sort_keys=True))
  else:
    print_summary(summary)
