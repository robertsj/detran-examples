# pyexamples/run_pin.py
#
# Solves a fully reflected C5G7 pin eigenproblem.
#
# Using 7x7 volume-conserving pin cell, DD SN, and a 8x8 QR quadrature,
# the reference kinf's are:
#   UO2       -- 1.325939322
#   MOX 4.3%  -- 1.184571821
#   MOX 7.0%  -- 1.135647780
#   MOX 8.7%  -- 1.159790118

import numpy as np
import time
from detran import *
try:
  from .pins_c5g7 import get_pins
  from .material_c5g7 import get_materials
  from .plot_utils import *
except ImportError:
  from pins_c5g7 import get_pins
  from material_c5g7 import get_materials
  from plot_utils import *

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
    ),
    eigen_monitor_level=2,
  )
  inp = eigenvalue_input(
    groups=7,
    dimension=2,
    equation="diffusion",
    boundary="reflect",
    inner=IterationSolver.source_iteration(
      max_iters=100,
      tolerance=1e-9,
      print_level=0,
      print_interval=10,
    ),
    outer=IterationSolver.gauss_seidel(
      max_iters=100,
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
      "eigen_solver_db": callow_db,
    },
  )
  #-----------------------------------------------------------------------------#
  # Material
  #-----------------------------------------------------------------------------#
  mat = get_materials()
  #-----------------------------------------------------------------------------#
  # Geometry
  #-----------------------------------------------------------------------------#
  pins = get_pins(7, True)
  mesh = pins[3].mesh()
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
    silo.initialize("c5g7pin.silo")
    silo.write_scalar_flux(state)
    silo.finalize()
  except :
    print("Silo error?")

if __name__ == "__main__":
  #Manager.initialize(sys.argv)
  run()
