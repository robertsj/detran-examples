# pyexamples/takeda/model.py
#
# Shared Takeda benchmark model data and solver builders.

import math
import time

import detran as dt


BOUNDARY = ("reflect", "vacuum", "reflect", "vacuum", "reflect", "vacuum")

CM_X = [0.0, 15.0, 20.0, 25.0]
CM_Y = [0.0, 5.0, 15.0, 25.0]
CM_Z = [0.0, 15.0, 25.0]

FM_X = [30, 10, 10]
FM_Y = [10, 20, 20]
FM_Z = [30, 20]

MATERIAL_MAP = [
  # Fueled region
  0, 2, 1,
  0, 1, 1,
  1, 1, 1,
  # Above fuel
  1, 2, 1,
  1, 1, 1,
  1, 1, 1,
]

CASES = {
  "inserted": {
    "name": "Model 1, inserted control",
    "case": 1,
    "reference_keff": 0.9624,
    "control": {
      "sigma_t": (8.52325e-02, 2.17460e-01),
      "sigma_s": (
        (6.77241e-02, 0.0),
        (6.45461e-05, 3.52358e-02),
      ),
    },
    "eigen_tolerance": 1.0e-4,
    "outer_tolerance": 1.0e-5,
  },
  "withdrawn": {
    "name": "Model 1, withdrawn control",
    "case": 2,
    "reference_keff": 0.9780,
    "control": {
      "sigma_t": (1.28407e-02, 1.20676e-02),
      "sigma_s": (
        (1.27700e-02, 0.0),
        (2.40997e-05, 1.07387e-02),
      ),
    },
    "eigen_tolerance": 1.0e-6,
    "outer_tolerance": 1.0e-6,
  },
}


def normalize_case(case):
  key = str(case).lower().replace("-", "_")
  aliases = {
    "1": "inserted",
    "case1": "inserted",
    "case_1": "inserted",
    "insert": "inserted",
    "inserted_control": "inserted",
    "2": "withdrawn",
    "case2": "withdrawn",
    "case_2": "withdrawn",
    "withdraw": "withdrawn",
    "void": "withdrawn",
    "withdrawn_control": "withdrawn",
  }
  key = aliases.get(key, key)
  if key not in CASES:
    raise ValueError("case must be 'inserted' or 'withdrawn'")
  return key


def _coarsen(divisions, mesh_divisor):
  mesh_divisor = int(mesh_divisor)
  if mesh_divisor <= 0:
    raise ValueError("mesh_divisor must be positive")
  return [max(1, int(math.ceil(value / mesh_divisor))) for value in divisions]


def build_mesh(mesh_divisor=1):
  """Build the 3-D Takeda Model 1 mesh.

  ``mesh_divisor=1`` gives the original fine mesh used by the historical
  example.  Larger values coarsen the fine divisions for quick MOC smoke runs.
  """

  return dt.Mesh3D(
    _coarsen(FM_X, mesh_divisor),
    _coarsen(FM_Y, mesh_divisor),
    _coarsen(FM_Z, mesh_divisor),
    CM_X,
    CM_Y,
    CM_Z,
    MATERIAL_MAP,
  )


def build_material(case="inserted"):
  case = normalize_case(case)
  control = CASES[case]["control"]

  material = dt.Material(3, 2, "Takeda Model 1")

  # 0. Core
  material.set_sigma_t(0, 0, 2.23775e-01)
  material.set_sigma_t(0, 1, 1.03864e+00)
  material.set_sigma_s(0, 0, 0, 1.92423e-01)
  material.set_sigma_s(0, 1, 0, 2.28253e-02)
  material.set_sigma_s(0, 1, 1, 8.80439e-01)
  material.set_sigma_f(0, 0, 9.09319e-03)
  material.set_sigma_f(0, 1, 2.90183e-01)
  material.set_chi(0, 0, 1.0)

  # 1. Reflector
  material.set_sigma_t(1, 0, 2.50367e-01)
  material.set_sigma_t(1, 1, 1.64482e+00)
  material.set_sigma_s(1, 0, 0, 1.93446e-01)
  material.set_sigma_s(1, 1, 0, 5.65042e-02)
  material.set_sigma_s(1, 1, 1, 1.62452e+00)

  # 2. Control rod or voided control channel
  material.set_sigma_t(2, 0, control["sigma_t"][0])
  material.set_sigma_t(2, 1, control["sigma_t"][1])
  material.set_sigma_s(2, 0, 0, control["sigma_s"][0][0])
  material.set_sigma_s(2, 0, 1, control["sigma_s"][0][1])
  material.set_sigma_s(2, 1, 0, control["sigma_s"][1][0])
  material.set_sigma_s(2, 1, 1, control["sigma_s"][1][1])

  material.finalize()
  return material


def build_legacy_input(case="inserted", equation="dd", quadrature=None):
  case = normalize_case(case)
  spec = CASES[case]
  solver_db = dt.SolverDatabase(
    linear=dt.LinearSolver.gmres(
      atol=1.0e-8,
      rtol=1.0e-8,
      max_iters=5000,
      restart=30,
      monitor_level=0,
    ),
    eigen_algorithm="gd",
    eigen_monitor_level=0,
  )
  return dt.eigenvalue_input(
    groups=2,
    dimension=3,
    equation=equation,
    boundary=BOUNDARY,
    inner=dt.IterationSolver.gmres(
      max_iters=1,
      tolerance=1.0e-6,
      print_level=0,
      print_interval=10,
      use_pc=False,
      pc_tolerance=1.0e-2,
      pc_ilu_levels=0,
      pc_type="none",
    ),
    outer=dt.IterationSolver.gmres(
      max_iters=0,
      tolerance=spec["outer_tolerance"],
      print_level=0,
      print_interval=2,
      upscatter_cutoff=0,
      use_pc=False,
    ),
    solver=dt.EigenSolver.arnoldi(
      max_iters=1000,
      tolerance=spec["eigen_tolerance"],
      monitor_level=2,
      print_interval=1,
    ),
    quadrature=quadrature or dt.Quadrature(polar=2, azimuth=2),
    solver_databases={
      "inner_solver_db": solver_db,
      "outer_solver_db": solver_db,
      "eigen_solver_db": solver_db,
    },
  )


def build_moc_problem(
  case="inserted",
  *,
  equation="scmoc",
  mesh_divisor=10,
  polar=1,
  azimuth=1,
  quadrature="u-u",
  spacing=5.0,
  inner_tolerance=1.0e-4,
  outer_tolerance=1.0e-4,
  eigen_tolerance=1.0e-4,
  inner_iters=40,
  outer_iters=40,
  eigen_iters=20,
  monitor_level=1,
  print_interval=1,
):
  mesh = build_mesh(mesh_divisor=mesh_divisor)
  material = build_material(case)
  problem = dt.EigenvalueProblem(
    mesh,
    material,
    equation=equation,
    boundary=BOUNDARY,
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
  return problem


def _try_write_input(prefix, inp, material, mesh):
  try:
    io = dt.IO_HDF5("{0}.h5".format(prefix))
    io.open()
    io.write(inp)
    io.write(material)
    io.write(mesh)
    io.close()
  except Exception:
    print("Error using HDF5---perhaps detran is not built with it.")


def _try_write_silo(prefix, mesh, state):
  try:
    silo = dt.SiloOutput(mesh)
    silo.initialize("{0}.silo".format(prefix))
    silo.write_mesh_map("MATERIAL")
    silo.write_scalar_flux(state)
    silo.finalize()
  except Exception:
    print("Error using Silo---perhaps detran is not built with it.")


def run_legacy(case="inserted", *, output_prefix=None, write_outputs=True):
  case = normalize_case(case)
  inp = build_legacy_input(case)
  mesh = build_mesh(mesh_divisor=1)
  material = build_material(case)
  prefix = output_prefix or "takeda_{0}".format(case)

  if write_outputs:
    _try_write_input(prefix, inp, material, mesh)

  start = time.time()
  solver = dt.Eigen3D(inp, material, mesh)
  solver.solve()
  elapsed = time.time() - start
  print("elapsed = ", elapsed)
  state = solver.state()

  if write_outputs:
    _try_write_silo(prefix, mesh, state)

  return {
    "case": case,
    "name": CASES[case]["name"],
    "keff": state.eigenvalue(),
    "reference_keff": CASES[case]["reference_keff"],
    "elapsed": elapsed,
    "state": state,
    "solver": solver,
  }
