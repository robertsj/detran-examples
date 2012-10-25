# pyexamples/run_core.py
#
# Solves the 2-D C5G7 benchmark problem.
#
# Using 7x7 volume-conserving pin cell, DD SN, and a QR quadrature
# of order 18, the reference keff is 1.186572165.

import numpy as np
import time
import sys
from detran import *
from core_c5g7 import get_core
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
inp.put_int("inner_max_iters",          1)
inp.put_dbl("inner_tolerance",          1e-3)
inp.put_int("inner_print_level",        0)
inp.put_int("inner_print_interval",     10)
#
inp.put_str("outer_solver",             "GS")
inp.put_int("outer_max_iters",          0)
inp.put_dbl("outer_tolerance",          1e-4)
inp.put_int("outer_print_level",        0)
inp.put_int("outer_print_interval",     1)
#
inp.put_str("eigen_solver",             "PI")
inp.put_int("eigen_max_iters",          10)
inp.put_dbl("eigen_tolerance",          1e-6)
inp.put_int("eigen_print_level",        2)
inp.put_int("eigen_print_interval",     1)
#
inp.put_str("bc_west",                  "reflect")
inp.put_str("bc_east",                  "vacuum")
inp.put_str("bc_south",                 "reflect")
inp.put_str("bc_north",                 "vacuum")
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
core = get_core(7, True)
mesh = core.mesh()
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
  silo.initialize("c5g7core.silo")
  silo.write_scalar_flux(state)
  silo.finalize()
except :
  print "Silo error?"
#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
Manager.finalize()

