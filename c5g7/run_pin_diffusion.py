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
inp = InputDB.Create()
#
inp.put_int("number_groups",            7)
#
inp.put_int("diffusion_max_iters",      1000)
inp.put_dbl("diffusion_tolerance",      1e-8)
#
inp.put_str("bc_left",                  "reflect")
inp.put_str("bc_right",                 "reflect")
inp.put_str("bc_bottom",                "reflect")
inp.put_str("bc_top",                   "reflect")
#-----------------------------------------------------------------------------#
# Material
#-----------------------------------------------------------------------------#
mat = get_materials()
for m in range(0, 7) :
  for g in range(0, 7) :
    mat.set_diff_coef(m, g, 1.0/(3.0*mat.sigma_t(m, g)))
mat.finalize()
#-----------------------------------------------------------------------------#
# Geometry
#-----------------------------------------------------------------------------#
pins = get_pins(7, True)
mesh = pins[0].mesh()
#-----------------------------------------------------------------------------#
# State
#-----------------------------------------------------------------------------#
state = State.Create(inp, mesh)
#-----------------------------------------------------------------------------#
# Execute
#-----------------------------------------------------------------------------#
Manager.initialize(sys.argv)
solver = DiffusionEigensolver(inp, mat, mesh, state)
t = time.time()
solver.solve()
print "elapsed = ", time.time()-t
#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
Manager.finalize()

