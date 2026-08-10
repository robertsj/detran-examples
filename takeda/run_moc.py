# pyexamples/takeda/run_moc.py
#
# Native 3-D MOC eigenvalue drivers for Takeda Model 1.

import argparse
import json
import math
import time

try:
  from .model import BOUNDARY, CASES, build_moc_problem, normalize_case
except ImportError:
  from model import BOUNDARY, CASES, build_moc_problem, normalize_case


def _material_counts(mesh):
  counts = {}
  for material_id in mesh.mesh_map("MATERIAL"):
    material_id = int(material_id)
    counts[material_id] = counts.get(material_id, 0) + 1
  return counts


def _flux_summary(result, mesh):
  groups = int(result.problem.groups)
  volumes = [float(mesh.volume(cell)) for cell in range(int(mesh.number_cells()))]
  total_volume = sum(volumes)
  summaries = []
  for group in range(groups):
    values = [float(value) for value in result.flux(group, numpy=False)]
    weighted = sum(value * volume for value, volume in zip(values, volumes))
    summaries.append({
      "group": group,
      "finite": all(math.isfinite(value) for value in values),
      "min": min(values),
      "max": max(values),
      "volume_average": weighted / total_volume if total_volume else 0.0,
    })
  return summaries


def run(**options):
  case = normalize_case(options.get("case", "inserted"))
  spec = CASES[case]

  start = time.time()
  problem = build_moc_problem(
    case,
    equation=options.get("equation", "scmoc"),
    mesh_divisor=int(options.get("mesh_divisor", 10)),
    polar=int(options.get("polar", 1)),
    azimuth=int(options.get("azimuth", 1)),
    quadrature=options.get("quadrature", "u-u"),
    spacing=float(options.get("spacing", 5.0)),
    inner_tolerance=float(options.get("inner_tolerance", 1.0e-4)),
    outer_tolerance=float(options.get("outer_tolerance", 1.0e-4)),
    eigen_tolerance=float(options.get("eigen_tolerance", 1.0e-4)),
    inner_iters=int(options.get("inner_iters", 40)),
    outer_iters=int(options.get("outer_iters", 40)),
    eigen_iters=int(options.get("eigen_iters", 20)),
    monitor_level=int(options.get("monitor_level", 1)),
    print_interval=int(options.get("print_interval", 1)),
  )
  result = problem.solve()
  elapsed = time.time() - start

  mesh = problem.mesh
  keff = float(result.keff)
  return {
    "case": case,
    "case_number": spec["case"],
    "name": spec["name"],
    "equation": problem.equation,
    "boundary": BOUNDARY,
    "reference_keff": spec["reference_keff"],
    "keff": keff,
    "delta_keff": keff - spec["reference_keff"],
    "finite": math.isfinite(keff),
    "mesh_divisor": int(options.get("mesh_divisor", 10)),
    "cells": int(mesh.number_cells()),
    "number_cells_x": int(mesh.number_cells_x()),
    "number_cells_y": int(mesh.number_cells_y()),
    "number_cells_z": int(mesh.number_cells_z()),
    "material_counts": _material_counts(mesh),
    "polar": int(options.get("polar", 1)),
    "azimuth": int(options.get("azimuth", 1)),
    "quadrature": options.get("quadrature", "u-u"),
    "spacing": float(options.get("spacing", 5.0)),
    "elapsed": elapsed,
    "sweeps": int(result.solver.number_sweeps())
      if hasattr(result.solver, "number_sweeps") else None,
    "flux": _flux_summary(result, mesh),
  }


def print_summary(summary):
  print("case: {0} ({1})".format(summary["case"], summary["name"]))
  print("equation: {0}".format(summary["equation"]))
  print("boundary(W,E,S,N,B,T): {0}".format(summary["boundary"]))
  print("mesh divisor: {0}".format(summary["mesh_divisor"]))
  print("mesh: {0} x {1} x {2} = {3} cells".format(
    summary["number_cells_x"],
    summary["number_cells_y"],
    summary["number_cells_z"],
    summary["cells"],
  ))
  print("material cell counts: {0}".format(summary["material_counts"]))
  print("quadrature: polar={0}, azimuth={1}, type={2}".format(
    summary["polar"],
    summary["azimuth"],
    summary["quadrature"],
  ))
  print("spacing: {0}".format(summary["spacing"]))
  print("reference keff: {:.8f}".format(summary["reference_keff"]))
  print("keff: {:.8f}".format(summary["keff"]))
  print("delta keff: {:+.8e}".format(summary["delta_keff"]))
  print("elapsed: {:.3f}s".format(summary["elapsed"]))
  print("sweeps: {0}".format(summary["sweeps"]))
  for flux in summary["flux"]:
    print(
      "group {0} flux finite/min/max/volume_average: {1} "
      "{2:.8e} {3:.8e} {4:.8e}".format(
        flux["group"],
        flux["finite"],
        flux["min"],
        flux["max"],
        flux["volume_average"],
      )
    )


def parse_args():
  parser = argparse.ArgumentParser(
    description="Run native 3-D MOC Takeda Model 1 eigenvalue cases.",
  )
  parser.add_argument("--case", choices=tuple(sorted(CASES)), default="inserted")
  parser.add_argument("--equation", choices=("scmoc", "ddmoc"), default="scmoc")
  parser.add_argument("--mesh-divisor", type=int, default=10,
                      help="coarsen original fine divisions by this factor")
  parser.add_argument("--reference-mesh", action="store_true",
                      help="use the original fine mesh, equivalent to --mesh-divisor 1")
  parser.add_argument("--polar", type=int, default=1)
  parser.add_argument("--azimuth", type=int, default=1)
  parser.add_argument("--quadrature", default="u-u")
  parser.add_argument("--spacing", type=float, default=5.0)
  parser.add_argument("--inner-tolerance", type=float, default=1.0e-4)
  parser.add_argument("--outer-tolerance", type=float, default=1.0e-4)
  parser.add_argument("--eigen-tolerance", type=float, default=1.0e-4)
  parser.add_argument("--inner-iters", type=int, default=40)
  parser.add_argument("--outer-iters", type=int, default=40)
  parser.add_argument("--eigen-iters", type=int, default=20)
  parser.add_argument("--monitor-level", type=int, default=1)
  parser.add_argument("--print-interval", type=int, default=1)
  parser.add_argument("--json", action="store_true")
  return parser.parse_args()


if __name__ == "__main__":
  args = parse_args()
  mesh_divisor = 1 if args.reference_mesh else args.mesh_divisor
  summary = run(
    case=args.case,
    equation=args.equation,
    mesh_divisor=mesh_divisor,
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
  )
  if args.json:
    print(json.dumps(summary, indent=2, sort_keys=True))
  else:
    print_summary(summary)
