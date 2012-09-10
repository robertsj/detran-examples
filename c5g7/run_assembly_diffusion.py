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
#
inp.put_int("number_groups",            7)
inp.put_int("diffusion_max_iters",      1000)
inp.put_dbl("diffusion_tolerance",      1e-6)
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
assemblies = get_assemblies(7, True)
mesh = assemblies[0].mesh()
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
# Plot
#-----------------------------------------------------------------------------#
from plot_utils import *
plot_flux(mesh, np.asarray(state.phi(0)))
#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
Manager.finalize()
