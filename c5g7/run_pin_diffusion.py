# pyexamples/run_pin.py
#
# Solves a fully reflected C5G7 pin eigenproblem using diffusion.
#

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
inp = InputDB()
inp.put_str("equation",                "diffusion")
inp.put_str("problem_type",                   "eigenvalue")
inp.put_int("number_groups",            7)
inp.put_str("bc_west",                        "reflect")
inp.put_str("bc_east",                        "reflect")
inp.put_str("bc_south",                       "reflect")
inp.put_str("bc_north",                       "reflect")
inp.put_int("eigen_max_iters",                1000)
inp.put_dbl("eigen_tolerance",                1e-9)

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
