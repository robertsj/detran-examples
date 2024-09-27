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
inp = InputDB()
inp.put_int("number_groups",                      7)
inp.put_int("dimension",                          2)
inp.put_str("equation",                           "diffusion")
inp.put_str("bc_west",                            "reflect")
inp.put_str("bc_east",                            "reflect")
inp.put_str("bc_south",                           "reflect")
inp.put_str("bc_north",                           "reflect")
inp.put_int("eigen_max_iters",                    1000)
## TODO: callow
inp.put_str("eigen_solver",                       "arnoldi")
db = InputDB("callow_db")
# outer gmres parameters
db.put_dbl("linear_solver_atol",                  1e-8);
db.put_dbl("linear_solver_rtol",                  1e-8);
db.put_str("linear_solver_type",                  "gmres");
db.put_int("linear_solver_maxit",                 5000);
db.put_int("linear_solver_gmres_restart",         30);
db.put_str("eigen_solver_type",                   "gd");
db.put_int("eigen_solver_monitor_level",          0);
db.put_int("linear_solver_monitor_level",         0);
inp.put_spdb("outer_solver_db", db)
inp.put_spdb("eigen_solver_db", db)
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
# Execute
#-----------------------------------------------------------------------------#
#Manager.initialize(sys.argv)
solver = Eigen2D(inp, mat, mesh)
t = time.time()
solver.solve()
print("elapsed = ", time.time()-t)
#-----------------------------------------------------------------------------#
# Plot
#-----------------------------------------------------------------------------#
from plot_utils import *
state = solver.state()
plot_flux(mesh, np.asarray(state.phi(0)))
#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
#Manager.finalize()
