# pyexamples/run_assembly.py
#
# Solves a fully reflected C5G7 assembly eigenproblem.
#
# Using 7x7 volume-conserving pin cell, DD SN, and a QR quadrature
# of order 50, the reference kinf's are:
#   UO2 -- 1.333964461
#   MOX -- 1.184571821

import numpy as np
import time
from detran import *
try:
  from .assemblies_c5g7 import get_assemblies
  from .pins_c5g7 import get_pins
  from .material_c5g7 import get_materials
except ImportError:
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
    boundary="reflect",
    inner=IterationSolver.source_iteration(
      max_iters=1,
      tolerance=1e-9,
      print_level=0,
      print_interval=10,
    ),
    outer=IterationSolver.gauss_seidel(
      max_iters=1,
      tolerance=1e-9,
      print_level=0,
      print_interval=1,
    ),
    solver=EigenSolver.power_iteration(
      max_iters=1000,
      tolerance=1e-9,
      monitor_level=2,
      print_interval=1,
      omega=1.0,
    ),
    quadrature=Quadrature(polar=5, azimuth=10),
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
  assemblies = get_assemblies(7, True)
  mesh = assemblies[1].mesh()
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
    silo.initialize("c5g7assembly.silo")
    silo.write_scalar_flux(state)
    silo.finalize()
  except :
    print("Silo error?")

if __name__ == "__main__":
  ##Manager.initialize(sys.argv)
  run()
  ##Manager.finalize()
