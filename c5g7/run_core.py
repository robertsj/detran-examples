# pyexamples/run_core.py
#
# Solves the 2-D C5G7 benchmark problem.
#
# Using 7x7 volume-conserving pin cell, DD SN, and a 3x6 QR quadrature,
# the reference keff is 1.186572165.
#
# By default, the problem is solved using generalized Davidson method
# using a multigroup, coarse-mesh, diffusion precondition (MGCMDSA).

import numpy as np
import time
from detran import *
try:
  from .core_c5g7 import get_core
  from .assemblies_c5g7 import get_assemblies
  from .pins_c5g7 import get_pins
  from .material_c5g7 import get_materials
except ImportError:
  from core_c5g7 import get_core
  from assemblies_c5g7 import get_assemblies
  from pins_c5g7 import get_pins
  from material_c5g7 import get_materials

def run() :

    callow_db = SolverDatabase(
      linear=LinearSolver.gmres(
        atol=1e-9,
        rtol=1e-9,
        max_iters=5000,
        restart=30,
        monitor_level=1,
        pc_type="ilu0",
      ),
      eigen_algorithm="gd",
      eigen_monitor_level=1,
    )
    inp = eigenvalue_input(
      groups=7,
      dimension=2,
      equation="dd",
      boundary=("reflect", "vacuum", "reflect", "vacuum"),
      inner=IterationSolver.gmres(
        max_iters=1,
        tolerance=1e-3,
        print_level=0,
        print_interval=10,
      ),
      outer=IterationSolver.gmres(
        max_iters=1000,
        tolerance=1e-8,
        print_level=0,
        print_interval=1,
        krylov_group_cutoff=0,
        pc_type="mgcmdsa",
      ),
      solver=EigenSolver(
        method="GD",
        algorithm="gd",
        max_iters=1000,
        tolerance=1e-6,
        monitor_level=2,
        print_interval=2,
      ),
      quadrature=Quadrature(type="asqr-asdr", polar=3, azimuth=6),
      solver_databases={
        "inner_solver_db": callow_db,
        "inner_pc_db": callow_db,
        "outer_solver_db": callow_db,
        "outer_pc_db": callow_db,
        "eigen_solver_db": callow_db,
        "eigen_solver_pc_db": callow_db,
      },
      input_overrides={
        "eigen_solver_pc_type": "mgcmdsa",
        "mgpc_coarse_mesh_level": 7,
        "mgpc_condensation_option": 1,
        "mgpc_cmdsa_use_smoothing": 1,
        "mgpc_cmdsa_smoothing_iters": 2,
        "mgpc_cmdsa_smoothing_relax": 0.8,
      },
    )

    mat = get_materials()
    core = get_core(7, True)
    mesh = core.mesh()

    start = time.time()
    solver = Eigen2D(inp, mat, mesh)
    solver.solve()
    print("elapsed = ", time.time() - start)


if __name__ == "__main__":
  #Manager.initialize(sys.argv)
  run()
  #Manager.finalize()
