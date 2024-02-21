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
import sys
from detran import *
from assemblies_c5g7 import get_assemblies
from pins_c5g7 import get_pins
from material_c5g7 import get_materials

def run() :

  #-----------------------------------------------------------------------------#
  # Input
  #-----------------------------------------------------------------------------#
  inp = InputDB()
  inp.put_str("equation",                       "dd")
  inp.put_str("problem_type",                   "eigenvalue")
  inp.put_int("number_groups",                  7)
  #
  inp.put_str("inner_solver",                   "SI")
  inp.put_int("inner_max_iters",                1)
  inp.put_dbl("inner_tolerance",                1e-9)
  inp.put_int("inner_print_level",              0)
  inp.put_int("inner_print_interval",           10)
  #
  inp.put_str("outer_solver",                   "GS")
  inp.put_int("outer_max_iters",                1)
  inp.put_dbl("outer_tolerance",                1e-9)
  inp.put_int("outer_print_level",              0)
  inp.put_int("outer_print_interval",           1)
  #
  inp.put_str("eigen_solver",                   "PI")
  inp.put_int("eigen_max_iters",                1000)
  inp.put_dbl("eigen_tolerance",                1e-9)
  inp.put_int("eigen_print_level",              2)
  inp.put_int("eigen_print_interval",           1)
  inp.put_dbl("eigen_pi_omega",                 1.0)
  #
  inp.put_str("bc_west",                        "reflect")
  inp.put_str("bc_east",                        "reflect")
  inp.put_str("bc_south",                       "reflect")
  inp.put_str("bc_north",                       "reflect")
  #
  inp.put_int("quad_number_polar_octant",       5)
  inp.put_int("quad_number_azimuth_octant",     10)
  #
  db = InputDB("callow_db")
  db.put_dbl("linear_solver_atol",              1e-9);
  db.put_dbl("linear_solver_rtol",              1e-9);
  db.put_str("linear_solver_type",              "petsc");
  db.put_int("linear_solver_maxit",             5000);
  db.put_int("linear_solver_gmres_restart",     30);
  db.put_int("linear_solver_monitor_level",     0);
  db.put_str("pc_type",                         "petsc_pc");
  db.put_str("petsc_pc_type",                   "lu");
  db.put_str("eigen_solver_type",               "slepc");
  db.put_int("eigen_solver_monitor_level",      2);
  inp.put_spdb("inner_solver_db",               db)
  inp.put_spdb("inner_pc_db",                   db)
  inp.put_spdb("outer_solver_db",               db)
  inp.put_spdb("eigen_solver_db",               db)
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
