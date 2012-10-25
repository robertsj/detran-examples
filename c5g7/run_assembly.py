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
from plot_utils import *
#-----------------------------------------------------------------------------#
# Initialize
#-----------------------------------------------------------------------------#
Manager.initialize(sys.argv)
#-----------------------------------------------------------------------------#
# Input
#-----------------------------------------------------------------------------#
inp = InputDB.Create()
inp.put_str("equation",                 "dd")
inp.put_str("problem_type",             "eigenvalue")
inp.put_int("number_groups",            7)
#
inp.put_str("inner_solver",             "SI")
inp.put_int("inner_max_iters",          10)
inp.put_dbl("inner_tolerance",          1e-3)
inp.put_int("inner_print_level",        0)
inp.put_int("inner_print_interval",     10)
#
inp.put_str("outer_solver",             "GS")
inp.put_int("outer_max_iters",          10)
inp.put_dbl("outer_tolerance",          1e-4)
inp.put_int("outer_print_level",        0)
inp.put_int("outer_print_interval",     1)
#
inp.put_str("eigen_solver",             "PI")
inp.put_int("eigen_max_iters",          100)
inp.put_dbl("eigen_tolerance",          1e-6)
inp.put_int("eigen_print_level",        2)
inp.put_int("eigen_print_interval",     1)
inp.put_dbl("eigen_pi_omega",           1.2)
#
inp.put_str("bc_west",                  "reflect")
inp.put_str("bc_east",                  "reflect")
inp.put_str("bc_south",                 "reflect")
inp.put_str("bc_north",                 "reflect")
#
inp.put_str("quad_type",                "quadruplerange")
inp.put_int("quad_order",               18)
#-----------------------------------------------------------------------------#
# Material
#-----------------------------------------------------------------------------#
mat = get_materials()
#-----------------------------------------------------------------------------#
# Geometry
#-----------------------------------------------------------------------------#
assemblies = get_assemblies(7, True)
mesh = assemblies[0].mesh()
#-----------------------------------------------------------------------------#
# Solve
#-----------------------------------------------------------------------------#
start = time.time()
solver = Eigen2D(inp, mat, mesh)
solver.solve()
print "elapsed = ", time.time() - start
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
  print "Silo error?"
#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
Manager.finalize()
