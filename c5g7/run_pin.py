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
# Input
#-----------------------------------------------------------------------------#
inp = InputDB.Create()
inp.put_str("equation",                 "dd")
inp.put_str("problem_type",             "eigenvalue")
#
inp.put_str("inner_solver",             "SI")
inp.put_int("inner_max_iters",          1000)
inp.put_dbl("inner_tolerance",          1e-8)
inp.put_int("inner_print_out",          0)
inp.put_int("inner_print_interval",     10)
#
inp.put_str("outer_solver",             "GS")
inp.put_int("outer_max_iters",          1000)
inp.put_dbl("outer_tolerance",          1e-8)
inp.put_int("outer_print_out",          0)
inp.put_int("outer_print_interval",     2)
# 
inp.put_str("eigen_solver",             "PI")
inp.put_int("eigen_max_iters",          1000)
inp.put_dbl("eigen_tolerance",          1e-8)
inp.put_int("eigen_print_out",          2)
inp.put_int("eigen_print_interval",     1)
#
inp.put_str("bc_left",                  "reflect")
inp.put_str("bc_right",                 "reflect")
inp.put_str("bc_bottom",                "reflect")
inp.put_str("bc_top",                   "reflect")
#
inp.put_str("quad_type",                "quadruplerange")
inp.put_int("quad_order",               2)
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
# Execute
#-----------------------------------------------------------------------------#
execute = Execute2D(sys.argv)
execute.initialize(inp, mat, mesh)
t = time.time()
execute.solve()
print "elapsed = ", time.time()-t
#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
execute.finalize()

