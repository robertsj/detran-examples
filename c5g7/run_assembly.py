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
# Input
#-----------------------------------------------------------------------------#
inp = InputDB.Create()
inp.put_str("equation",                 "dd")
inp.put_str("problem_type",             "eigenvalue")
#
inp.put_str("inner_solver",             "SI")
inp.put_int("inner_max_iters",          1)
inp.put_dbl("inner_tolerance",          1e-8)
inp.put_int("inner_print_out",          0)
inp.put_int("inner_print_interval",     10)
inp.put_int("inner_use_pc",             1)
#
inp.put_str("outer_solver",             "GS")
inp.put_int("outer_max_iters",          0)
inp.put_dbl("outer_tolerance",          1e-8)
inp.put_int("outer_print_out",          1)
inp.put_int("outer_print_interval",     2)
inp.put_int("outer_use_pc",             0)
inp.put_dbl("outer_pc_tolerance",       1e-3)
inp.put_int("outer_upscatter_cutoff",   0)
# 
inp.put_str("eigen_solver",             "PI")
inp.put_int("eigen_max_iters",          100)
inp.put_dbl("eigen_tolerance",          1e-8)
inp.put_int("eigen_print_out",          2)
inp.put_int("eigen_print_interval",     1)
#
#inp.put_str("bc_left",                  "reflect")
#inp.put_str("bc_right",                 "reflect")
#inp.put_str("bc_bottom",                "reflect")
#inp.put_str("bc_top",                   "reflect")
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
assemblies = get_assemblies(7, True)
mesh = assemblies[0].mesh()
#-----------------------------------------------------------------------------#
# Execute
#-----------------------------------------------------------------------------#
execute = Execute2D(sys.argv)
execute.initialize(inp, mat, mesh)
t = time.time()
execute.solve()
print "elapsed = ", time.time()-t
#-----------------------------------------------------------------------------#
# Plot
#-----------------------------------------------------------------------------#
#from plot_utils import *
#state = execute.get_state()
#plot_flux(mesh, np.asarray(state.phi(0)))
#try:
#  silo = SiloOutput(inp, mesh)
#  silo.initialize()
#  silo.write_flux(state)
#  silo.finalize()
#except:
#  print "no silo"
#rates = ReactionRates(mat, mesh, state)
#pinpower = rates.region_power("PINS") 
#print np.asarray(pinpower)
#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
execute.finalize()

