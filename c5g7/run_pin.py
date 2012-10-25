# pyexamples/run_pin.py
#
# Solves a fully reflected C5G7 pin eigenproblem.
#
# Using 7x7 volume-conserving pin cell, DD SN, and a QR quadrature
# of order 50, the reference kinf's are:
#   UO2       -- 1.325902846
#   MOX 4.3%  -- 1.184571821
#   MOC 7.0%  -- 1.135647780
#   MOX 8.7%  -- 1.159790118

import numpy as np
import time
import sys
from detran import *
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
inp.put_int("eigen_max_iters",          170)
inp.put_dbl("eigen_tolerance",          1e-9)
inp.put_int("eigen_print_level",        2)
inp.put_int("eigen_print_interval",     1)
inp.put_dbl("eigen_pi_omega",           1.0)
#
inp.put_str("bc_west",                  "reflect")
inp.put_str("bc_east",                  "reflect")
inp.put_str("bc_south",                 "reflect")
inp.put_str("bc_north",                 "reflect")
#
inp.put_str("quad_type",                "quadruplerange")
inp.put_int("quad_order",               50)
#-----------------------------------------------------------------------------#
# Material
#-----------------------------------------------------------------------------#
mat = get_materials()
#-----------------------------------------------------------------------------#
# Geometry
#-----------------------------------------------------------------------------#
pins = get_pins(7, True)
mesh = pins[0].mesh()
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
  silo.initialize("c5g7pin.silo")
  silo.write_scalar_flux(state)
  silo.finalize()
except :
  print "Silo error?"
#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
Manager.finalize()
