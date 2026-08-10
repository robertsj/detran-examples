# pyexamples/c5g7/moc_driver.py
#
# Shared native 2-D MOC driver support for the C5G7 pin, assembly, and core
# examples.

import argparse
import math
import time

import detran as dt

try:
  from .assemblies_c5g7 import get_assemblies
  from .core_c5g7 import get_core
  from .material_c5g7 import get_materials
  from .pins_c5g7 import get_pins
  from .plot_utils import write_moc_flux_svgs
  from .reference import assembly_eigenvalues, eigenvalue, pin_eigenvalues
except ImportError:
  from assemblies_c5g7 import get_assemblies
  from core_c5g7 import get_core
  from material_c5g7 import get_materials
  from pins_c5g7 import get_pins
  from plot_utils import write_moc_flux_svgs
  from reference import assembly_eigenvalues, eigenvalue, pin_eigenvalues


MOC_OPTION_NAMES = (
  "index",
  "pin_mesh",
  "polar",
  "azimuth",
  "quadrature",
  "spacing",
  "equation",
  "inner_tolerance",
  "outer_tolerance",
  "eigen_tolerance",
  "inner_iters",
  "outer_iters",
  "eigen_iters",
  "monitor_level",
  "print_interval",
  "plot_dir",
  "plot_groups",
  "plot_style",
  "plot_width",
  "plot_palette",
)

MOC_PROBLEM_OPTION_NAMES = (
  "polar",
  "azimuth",
  "quadrature",
  "spacing",
  "equation",
  "inner_tolerance",
  "outer_tolerance",
  "eigen_tolerance",
  "inner_iters",
  "outer_iters",
  "eigen_iters",
  "monitor_level",
  "print_interval",
)

DEFAULT_MOC_OPTIONS = {
  "index": 0,
  "pin_mesh": 7,
  "polar": 1,
  "azimuth": 4,
  "quadrature": "u-u",
  "spacing": 0.25,
  "equation": "scmoc",
  "inner_tolerance": 1.0e-5,
  "outer_tolerance": 1.0e-5,
  "eigen_tolerance": 1.0e-4,
  "inner_iters": 80,
  "outer_iters": 80,
  "eigen_iters": 60,
  "monitor_level": 1,
  "print_interval": 10,
  "plot_dir": "moc_plots",
  "plot_groups": None,
  "plot_style": "tracks",
  "plot_width": 900,
  "plot_palette": "viridis",
}

CASE_MOC_OPTIONS = {
  "core": {
    "plot_groups": (0,),
    "plot_style": "cells",
  },
}

C5G7_CORE_BOUNDARY = ("reflect", "vacuum", "reflect", "vacuum")

SUMMARY_KEYS = (
  "case",
  "index",
  "cells",
  "width_x",
  "width_y",
  "boundary",
  "polar",
  "azimuth",
  "spacing",
  "equation",
  "keff",
  "finite",
  "reference_name",
  "reference_keff",
  "delta_keff",
  "elapsed",
  "sweeps",
  "plot_files",
)


def parse_plot_groups(spec):
  """Parse a comma-separated, zero-based group selection."""
  if spec is None:
    return None
  spec = str(spec).strip()
  if not spec or spec.lower() == "all":
    return "all"

  groups = []
  for part in spec.split(","):
    part = part.strip()
    if not part:
      continue
    if "-" in part:
      left, right = part.split("-", 1)
      start = int(left)
      stop = int(right)
      if stop < start:
        raise argparse.ArgumentTypeError("plot group ranges must be increasing")
      groups.extend(range(start, stop + 1))
    else:
      groups.append(int(part))

  if not groups:
    raise argparse.ArgumentTypeError("plot groups must not be empty")
  if any(group < 0 for group in groups):
    raise argparse.ArgumentTypeError("plot groups must be nonnegative")

  return tuple(dict.fromkeys(groups))


def get_moc_case(case, index=0, pin_mesh=7):
  """Return the mesh and reference data for a C5G7 MOC case."""
  if case == "pin":
    pins = get_pins(pin_mesh, True)
    if index < 0 or index >= len(pins):
      raise ValueError("pin index must be between 0 and 6")
    return pins[index].mesh(), pin_eigenvalues().get(index)

  if case == "assembly":
    assemblies = get_assemblies(pin_mesh, True)
    if index < 0 or index >= len(assemblies):
      raise ValueError("assembly index must be between 0 and 2")
    return assemblies[index].mesh(), assembly_eigenvalues().get(index)

  if case == "core":
    if index not in (None, 0):
      raise ValueError("core index must be 0")
    return get_core(pin_mesh, True).mesh(), ("C5G7 core", eigenvalue())

  raise ValueError("case must be 'pin', 'assembly', or 'core'")


def get_moc_boundary(case):
  """Return the boundary condition for a C5G7 MOC case."""
  if case in ("pin", "assembly"):
    return "reflect"
  if case == "core":
    return C5G7_CORE_BOUNDARY
  raise ValueError("case must be 'pin', 'assembly', or 'core'")


def build_moc_problem(
  mesh,
  material,
  *,
  polar=1,
  azimuth=4,
  quadrature="u-u",
  spacing=0.25,
  equation="scmoc",
  inner_tolerance=1.0e-5,
  outer_tolerance=1.0e-5,
  eigen_tolerance=1.0e-4,
  inner_iters=80,
  outer_iters=80,
  eigen_iters=60,
  monitor_level=1,
  print_interval=10,
  boundary="reflect",
):
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
      "tracker_maximum_spacing": spacing,
      "tracker_minimum_number": 1,
      "tracker_spatial_quad_type": "uniform",
      "inner_print_level": 0,
      "outer_print_level": 0,
    },
  )


def summarize_moc_result(case, index, mesh, result, reference, options, elapsed):
  summary = {
    "case": case,
    "index": index,
    "cells": mesh.number_cells(),
    "width_x": mesh.total_width_x(),
    "width_y": mesh.total_width_y(),
    "boundary": options["boundary"],
    "polar": options["polar"],
    "azimuth": options["azimuth"],
    "spacing": options["spacing"],
    "equation": options["equation"],
    "keff": result.keff,
    "finite": math.isfinite(result.keff),
    "elapsed": elapsed,
  }
  if hasattr(result.solver, "number_sweeps"):
    summary["sweeps"] = result.solver.number_sweeps()
  if reference is not None:
    summary["reference_name"] = reference[0]
    summary["reference_keff"] = reference[1]
    summary["delta_keff"] = result.keff - reference[1]
  return summary


def _plot_groups(options, material):
  if options["plot_groups"] == "all":
    groups = tuple(range(int(material.number_groups())))
  elif options["plot_groups"] is not None:
    groups = tuple(int(group) for group in options["plot_groups"])
  else:
    groups = tuple(range(int(material.number_groups())))
  invalid = [group for group in groups if group >= int(material.number_groups())]
  if invalid:
    raise ValueError(
      "plot group ids must be less than {0}: {1}".format(
        material.number_groups(),
        invalid,
      )
    )
  return groups


def run_moc_case(case, **options):
  defaults = dict(DEFAULT_MOC_OPTIONS)
  defaults.update(CASE_MOC_OPTIONS.get(case, {}))
  defaults["boundary"] = get_moc_boundary(case)
  options = dict(defaults, **options)
  material = get_materials()
  index = options["index"]
  mesh, reference = get_moc_case(
    case,
    index=index,
    pin_mesh=options["pin_mesh"],
  )
  problem_options = {name: options[name] for name in MOC_PROBLEM_OPTION_NAMES}
  problem_options["boundary"] = options["boundary"]
  problem = build_moc_problem(mesh, material, **problem_options)

  start = time.time()
  result = problem.solve()
  elapsed = time.time() - start

  summary = summarize_moc_result(case, index, mesh, result, reference, options, elapsed)
  if options["plot_dir"]:
    summary["plot_files"] = write_moc_flux_svgs(
      mesh,
      result,
      options["plot_dir"],
      "c5g7_{0}_{1}_{2}".format(case, index, options["equation"]),
      groups=_plot_groups(options, material),
      width=options["plot_width"],
      palette=options["plot_palette"],
      style=options["plot_style"],
    )
  return summary


def print_summary(summary):
  for key in SUMMARY_KEYS:
    if key in summary:
      print(f"{key}: {summary[key]}")


def add_moc_arguments(
  parser,
  *,
  include_case=False,
  default_case="assembly",
  default_index=0,
):
  if include_case:
    parser.add_argument("--case", choices=("pin", "assembly", "core"),
                        default=default_case)
  parser.add_argument("--index", type=int, default=default_index,
                      help="pin or assembly index; core requires 0")
  parser.add_argument("--pin-mesh", type=int, default=7,
                      help="pin-cell mesh subdivisions per dimension")
  parser.add_argument("--polar", type=int, default=1,
                      help="product polar angles per octant")
  parser.add_argument("--azimuth", type=int, default=4,
                      help="product azimuthal angles per octant")
  parser.add_argument("--quadrature", default="u-u",
                      help="product quadrature type")
  parser.add_argument("--spacing", type=float, default=0.25,
                      help="maximum average track spacing")
  parser.add_argument("--equation", choices=("scmoc", "ddmoc"), default="scmoc",
                      help="native MOC spatial discretization")
  parser.add_argument("--inner-tolerance", type=float, default=1.0e-5)
  parser.add_argument("--outer-tolerance", type=float, default=1.0e-5)
  parser.add_argument("--eigen-tolerance", type=float, default=1.0e-4)
  parser.add_argument("--inner-iters", type=int, default=80)
  parser.add_argument("--outer-iters", type=int, default=80)
  parser.add_argument("--eigen-iters", type=int, default=60)
  parser.add_argument("--monitor-level", type=int, default=1)
  parser.add_argument("--print-interval", type=int, default=10)
  parser.add_argument("--plot-dir", default="moc_plots",
                      help="directory for scalar-flux SVG plots")
  parser.add_argument("--no-plots", action="store_true",
                      help="skip SVG flux plot generation")
  parser.add_argument("--plot-groups", type=parse_plot_groups, default=None,
                      metavar="GROUPS",
                      help="zero-based group list/ranges, e.g. 0,2-4 or all")
  parser.add_argument("--plot-style", choices=("tracks", "cells"), default=None,
                      help="render flux as track segments or mesh cells")
  parser.add_argument("--plot-width", type=int, default=900,
                      help="SVG width in pixels")
  parser.add_argument("--plot-palette", default="viridis",
                      help="continuous SVG palette or matplotlib colormap name")
  return parser


def parse_moc_args(
  description,
  *,
  include_case=False,
  default_case="assembly",
  default_index=0,
):
  parser = argparse.ArgumentParser(description=description)
  add_moc_arguments(
    parser,
    include_case=include_case,
    default_case=default_case,
    default_index=default_index,
  )
  return parser.parse_args()


def moc_options_from_args(args):
  options = {}
  for name in MOC_OPTION_NAMES:
    value = getattr(args, name)
    if value is not None:
      options[name] = value
  if getattr(args, "no_plots", False):
    options["plot_dir"] = None
  return options


def run_from_args(args, *, case=None):
  selected_case = case or args.case
  return run_moc_case(selected_case, **moc_options_from_args(args))
