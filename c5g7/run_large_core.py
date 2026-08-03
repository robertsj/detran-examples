# pyexamples/run_core.py
#
# Solves the 2-D C5G7 benchmark problem.
#
# Using 7x7 volume-conserving pin cell, DD SN, and a 3x6 QR quadrature,
# the reference keff is 1.186572165.

import numpy as np
import time
from detran import *
try:
  from .core_c5g7 import get_large_core
  from .assemblies_c5g7 import get_assemblies
  from .pins_c5g7 import get_pins
  from .material_c5g7 import get_materials
except ImportError:
  from core_c5g7 import get_large_core
  from assemblies_c5g7 import get_assemblies
  from pins_c5g7 import get_pins
  from material_c5g7 import get_materials

def run() :
  #-----------------------------------------------------------------------------#
  # Input 
  #-----------------------------------------------------------------------------#
  callow_db = SolverDatabase(
    linear=LinearSolver.gmres(
      atol=1e-9,
      rtol=1e-9,
      max_iters=5000,
      restart=30,
      monitor_level=0,
      pc_type="ilu0",
      pc_side=1,
    ),
    eigen_algorithm="gd",
    eigen_monitor_level=2,
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
      max_iters=0,
      tolerance=1e-4,
      print_level=0,
      print_interval=1,
    ),
    solver=EigenSolver.power_iteration(
      max_iters=1000,
      tolerance=1e-6,
      monitor_level=2,
      print_interval=1,
    ),
    quadrature=Quadrature(type="quadruplerange", polar=3, azimuth=6),
    solver_databases={
      "inner_solver_db": callow_db,
      "inner_pc_db": callow_db,
      "outer_solver_db": callow_db,
    },
  )
  #-----------------------------------------------------------------------------#
  # Material
  #-----------------------------------------------------------------------------#
  mat = get_materials()
  #-----------------------------------------------------------------------------#
  # Geometry
  #-----------------------------------------------------------------------------#
  core = get_large_core(7, True)
  mesh = core.mesh()

  #-----------------------------------------------------------------------------#
  # Solve
  #-----------------------------------------------------------------------------#
  start = time.time()
  solver = Eigen2D(inp, mat, mesh)
  solver.solve()
  print("elapsed = ", time.time() - start)
  #-----------------------------------------------------------------------------#
  # Plot
  #-----------------------------------------------------------------------------#
  try :
    state = solver.state()
    silo = SiloOutput(mesh)
    silo.initialize("c5g7largecore.silo")
    silo.write_scalar_flux(state)
    silo.write_mesh_map("MATERIAL")
    silo.finalize()
  except :
    print("Silo error?")

if __name__ == "__main__":
  #Manager.initialize(sys.argv)
  run()
  #Manager.finalize()
