# pyexamples/kobayashi/run_moc.py
#
# Native 3-D MOC fixed-source drivers for the Kobayashi void benchmarks.

import argparse
import json
import math
import time

import detran as dt


BOUNDARY = ("reflect", "vacuum", "reflect", "vacuum", "reflect", "vacuum")


PROBLEMS = {
  1: {
    "name": "volume source in cubic shield",
    "widths": {
      "x": (10.0, 40.0, 50.0),
      "y": (10.0, 40.0, 50.0),
      "z": (10.0, 40.0, 50.0),
    },
    "divisions": {
      "x": (1, 4, 5),
      "y": (1, 4, 5),
      "z": (1, 4, 5),
    },
    "materials": (
      (
        (2, 0, 1),
        (0, 0, 1),
        (1, 1, 1),
      ),
      (
        (0, 0, 1),
        (0, 0, 1),
        (1, 1, 1),
      ),
      (
        (1, 1, 1),
        (1, 1, 1),
        (1, 1, 1),
      ),
    ),
    "lines": {
      "A": ((5.0, 5.0, 5.0), (5.0, 95.0, 5.0), 10),
      "B": ((5.0, 5.0, 5.0), (95.0, 95.0, 95.0), 10),
      "C": ((5.0, 55.0, 5.0), (95.0, 55.0, 5.0), 10),
    },
  },
  2: {
    "name": "straight duct",
    "widths": {
      "x": (10.0, 50.0),
      "y": (10.0, 90.0),
      "z": (10.0, 50.0),
    },
    "divisions": {
      "x": (1, 5),
      "y": (1, 9),
      "z": (1, 5),
    },
    "materials": (
      (
        (2, 1),
        (0, 1),
      ),
      (
        (1, 1),
        (1, 1),
      ),
    ),
    "lines": {
      "A": ((5.0, 5.0, 5.0), (5.0, 95.0, 5.0), 10),
      "B": ((5.0, 95.0, 5.0), (55.0, 95.0, 5.0), 6),
    },
  },
  3: {
    "name": "dog-legged duct",
    "widths": {
      "x": (10.0, 20.0, 10.0, 20.0),
      "y": (10.0, 40.0, 10.0, 40.0),
      "z": (10.0, 20.0, 10.0, 20.0),
    },
    "divisions": {
      "x": (1, 2, 1, 2),
      "y": (1, 4, 1, 4),
      "z": (1, 2, 1, 2),
    },
    "materials": (
      (
        (2, 1, 1, 1),
        (0, 1, 1, 1),
        (0, 0, 0, 1),
        (1, 1, 1, 1),
      ),
      (
        (1, 1, 1, 1),
        (1, 1, 1, 1),
        (1, 1, 0, 1),
        (1, 1, 1, 1),
      ),
      (
        (1, 1, 1, 1),
        (1, 1, 1, 1),
        (1, 1, 0, 1),
        (1, 1, 0, 1),
      ),
      (
        (1, 1, 1, 1),
        (1, 1, 1, 1),
        (1, 1, 1, 1),
        (1, 1, 1, 1),
      ),
    ),
    "lines": {
      "A": ((5.0, 5.0, 5.0), (5.0, 95.0, 5.0), 10),
      "B": ((5.0, 55.0, 5.0), (55.0, 55.0, 5.0), 6),
      "C": ((5.0, 95.0, 35.0), (55.0, 95.0, 35.0), 6),
    },
  },
}


def _edges(widths):
  total = 0.0
  edges = [total]
  for width in widths:
    total += float(width)
    edges.append(total)
  return edges


def _refined(divisions, refine):
  factor = 2 ** int(refine)
  return [int(value) * factor for value in divisions]


def _material_counts(mesh):
  counts = {}
  for material_id in mesh.mesh_map("MATERIAL"):
    material_id = int(material_id)
    counts[material_id] = counts.get(material_id, 0) + 1
  return counts


def build_material():
  material = dt.Material(3, 1, "Kobayashi")
  for material_id, sigma_t, sigma_s in (
    (0, 1.0e-4, 0.5e-4),
    (1, 0.1, 0.05),
    (2, 0.1, 0.05),
  ):
    material.set_sigma_t(material_id, 0, sigma_t)
    material.set_sigma_s(material_id, 0, 0, sigma_s)
  material.compute_sigma_a()
  material.finalize()
  return material


def build_mesh(spec, refine=0):
  widths = spec["widths"]
  divisions = spec["divisions"]
  return dt.CartesianMesh.from_coarse(
    x=_edges(widths["x"]),
    xdiv=_refined(divisions["x"], refine),
    y=_edges(widths["y"]),
    ydiv=_refined(divisions["y"], refine),
    z=_edges(widths["z"]),
    zdiv=_refined(divisions["z"], refine),
    materials=spec["materials"],
  )


def build_source(mesh, quadrature):
  source_map = [
    1 if int(material_id) == 2 else 0
    for material_id in mesh.mesh_map("MATERIAL")
  ]
  spectra = dt.vec2_dbl([[0.0], [1.0]])
  return dt.IsotropicSource(
    1,
    mesh,
    spectra,
    dt.vec_int(source_map),
    quadrature,
  )


def build_problem(
  mesh,
  material,
  *,
  equation="scmoc",
  polar=1,
  azimuth=1,
  quadrature="u-u",
  spacing=10.0,
  tolerance=1.0e-6,
  max_iters=80,
):
  return dt.FixedSourceProblem(
    mesh,
    material,
    equation=equation,
    boundary=BOUNDARY,
    inner=dt.IterationSolver.source_iteration(
      tolerance=tolerance,
      max_iters=max_iters,
      print_level=0,
    ),
    outer=dt.IterationSolver.gauss_seidel(
      tolerance=tolerance,
      max_iters=max_iters,
      print_level=0,
    ),
    quadrature=dt.Quadrature.product(
      polar=polar,
      azimuth=azimuth,
      type=quadrature,
    ),
    input_overrides={
      "tracker_maximum_spacing": spacing,
      "tracker_minimum_number": 1,
      "tracker_spatial_quad_type": "uniform",
      "inner_print_level": 0,
      "outer_print_level": 0,
    },
  )


def _flux_summary(values, volumes):
  total_volume = sum(volumes)
  finite = [value for value in values if math.isfinite(value)]
  weighted = sum(value * volume for value, volume in zip(values, volumes))
  return {
    "finite": len(finite) == len(values),
    "min": min(values),
    "max": max(values),
    "volume_average": weighted / total_volume if total_volume else 0.0,
  }


def run(**options):
  case = int(options.get("case", 1))
  spec = PROBLEMS[case]
  material = build_material()
  mesh = build_mesh(spec, refine=int(options.get("refine", 0)))
  problem = build_problem(
    mesh,
    material,
    equation=options.get("equation", "scmoc"),
    polar=int(options.get("polar", 1)),
    azimuth=int(options.get("azimuth", 1)),
    quadrature=options.get("quadrature", "u-u"),
    spacing=float(options.get("spacing", 10.0)),
    tolerance=float(options.get("tolerance", 1.0e-6)),
    max_iters=int(options.get("max_iters", 80)),
  )

  start = time.time()
  solver = problem.prepare_raw_solver(
    setup=lambda raw_solver: raw_solver.set_source(
      build_source(mesh, raw_solver.quadrature())
    )
  )
  ok = solver.solve()
  if ok is False:
    raise RuntimeError("Kobayashi fixed-source solve failed")
  elapsed = time.time() - start

  phi = [float(value) for value in solver.state().phi(0)]
  volumes = [float(mesh.volume(cell)) for cell in range(int(mesh.number_cells()))]
  return {
    "case": case,
    "name": spec["name"],
    "equation": options.get("equation", "scmoc"),
    "boundary": BOUNDARY,
    "cells": int(mesh.number_cells()),
    "number_cells_x": int(mesh.number_cells_x()),
    "number_cells_y": int(mesh.number_cells_y()),
    "number_cells_z": int(mesh.number_cells_z()),
    "material_counts": _material_counts(mesh),
    "source_cells": _material_counts(mesh).get(2, 0),
    "polar": int(options.get("polar", 1)),
    "azimuth": int(options.get("azimuth", 1)),
    "spacing": float(options.get("spacing", 10.0)),
    "elapsed": elapsed,
    "sweeps": int(solver.number_sweeps()) if hasattr(solver, "number_sweeps") else None,
    "flux": _flux_summary(phi, volumes),
    "sample_lines": spec["lines"],
  }


def print_summary(summary):
  print("case: {0} ({1})".format(summary["case"], summary["name"]))
  print("equation: {0}".format(summary["equation"]))
  print("boundary(W,E,S,N,B,T): {0}".format(summary["boundary"]))
  print("mesh: {0} x {1} x {2} = {3} cells".format(
    summary["number_cells_x"],
    summary["number_cells_y"],
    summary["number_cells_z"],
    summary["cells"],
  ))
  print("material cell counts: {0}".format(summary["material_counts"]))
  print("source cells: {0}".format(summary["source_cells"]))
  print("quadrature: polar={0}, azimuth={1}".format(
    summary["polar"],
    summary["azimuth"],
  ))
  print("spacing: {0}".format(summary["spacing"]))
  print("elapsed: {:.3f}s".format(summary["elapsed"]))
  print("sweeps: {0}".format(summary["sweeps"]))
  print("flux finite: {0}".format(summary["flux"]["finite"]))
  print("flux min/max/volume_average: {:.8e} {:.8e} {:.8e}".format(
    summary["flux"]["min"],
    summary["flux"]["max"],
    summary["flux"]["volume_average"],
  ))
  print("sample lines: {0}".format(tuple(sorted(summary["sample_lines"]))))


def parse_args():
  parser = argparse.ArgumentParser(
    description="Run a native 3-D MOC Kobayashi fixed-source benchmark smoke.",
  )
  parser.add_argument("--case", type=int, choices=tuple(sorted(PROBLEMS)), default=1)
  parser.add_argument("--equation", choices=("scmoc", "ddmoc"), default="scmoc")
  parser.add_argument("--refine", type=int, default=0,
                      help="uniformly multiply coarse interval divisions by 2**refine")
  parser.add_argument("--polar", type=int, default=1)
  parser.add_argument("--azimuth", type=int, default=1)
  parser.add_argument("--quadrature", default="u-u")
  parser.add_argument("--spacing", type=float, default=10.0)
  parser.add_argument("--tolerance", type=float, default=1.0e-6)
  parser.add_argument("--max-iters", type=int, default=80)
  parser.add_argument("--json", action="store_true")
  return parser.parse_args()


if __name__ == "__main__":
  args = parse_args()
  summary = run(
    case=args.case,
    equation=args.equation,
    refine=args.refine,
    polar=args.polar,
    azimuth=args.azimuth,
    quadrature=args.quadrature,
    spacing=args.spacing,
    tolerance=args.tolerance,
    max_iters=args.max_iters,
  )
  if args.json:
    print(json.dumps(summary, indent=2, sort_keys=True))
  else:
    print_summary(summary)
