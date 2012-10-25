# pyexamples/takeda/takeda_withdrawn.py
#
# Takeda benchmark, inserted control.
# Reference keff ~ 0.9780

import numpy as np
import time
import sys
from detran import *

# Initialize
Manager.initialize(sys.argv)

# Input
inp = InputDB.Create()
inp.put_str("problem_type",             "eigenvalue")
inp.put_int("number_groups",            2)
inp.put_int("dimension",                3)
inp.put_str("equation",                 "dd")
# Inner solver
inp.put_str("inner_solver",             "SI")
inp.put_int("inner_max_iters",          1)
inp.put_dbl("inner_tolerance",          1e-6)
inp.put_int("inner_print_out",          0)
inp.put_int("inner_print_interval",     10)
inp.put_int("inner_use_pc",             0)
inp.put_dbl("inner_pc_tolerance",       1e-2)
inp.put_int("inner_pc_ilu_levels",      0)
inp.put_str("inner_pc_type",            "none")
# Outer solver
inp.put_str("outer_solver",             "GS")
inp.put_int("outer_max_iters",          0)
inp.put_dbl("outer_tolerance",          1e-5)
inp.put_int("outer_print_out",          0)
inp.put_int("outer_print_interval",     2)
inp.put_int("outer_upscatter_cutoff",   0) 
inp.put_int("outer_use_pc",             0)
# Eigensolver
inp.put_str("eigen_solver",             "SLEPc")
inp.put_int("eigen_max_iters",          1000)
inp.put_dbl("eigen_tolerance",          1e-5)
inp.put_int("eigen_print_out",          2)
inp.put_int("eigen_print_interval",     1)
# Boundary
inp.put_str("bc_west",                  "reflect")
inp.put_str("bc_east",                  "vacuum")
inp.put_str("bc_south",                 "reflect")
inp.put_str("bc_north",                 "vacuum")
inp.put_str("bc_bottom",                "reflect")
inp.put_str("bc_top",                   "vacuum")
#
inp.put_str("quad_type",                "quadruplerange")
inp.put_int("quad_order",               50)

# Mesh
cm_x = [0.0, 15.0, 20.0, 25.0]
cm_y = [0.0,  5.0, 15.0, 25.0]
cm_z = [0.0, 15.0, 25.0]
fm_x = [30, 10, 10]
fm_y = [10, 20, 20]
fm_z = [30, 20]
mt   = [# Fueled region
        0, 2, 1, 
        0, 1, 1,
        1, 1, 1, 
        # Above fuel
        1, 2, 1, 
        1, 1, 1,
        1, 1, 1]

mesh = Mesh3D.Create(fm_x, fm_y, fm_z, cm_x, cm_y, cm_z, mt)

# Material
mat = Material.Create(3, 2, False)
# 0. Core
mat.set_sigma_t(0, 0,    2.23775e-01)
mat.set_sigma_t(0, 1,    1.03864e-00)
mat.set_sigma_s(0, 0, 0, 1.92423e-01)
mat.set_sigma_s(0, 1, 0, 2.28253e-02)
mat.set_sigma_s(0, 1, 1, 8.80439e-01)
mat.set_sigma_f(0, 0,    9.09319e-03)
mat.set_sigma_f(0, 1,    2.90183e-01)
mat.set_chi(0, 0,        1.0)
# 1. Reflector
mat.set_sigma_t(1, 0,    2.50367e-01)
mat.set_sigma_t(1, 1,    1.64482e-00)
mat.set_sigma_s(1, 0, 0, 1.93446e-01)
mat.set_sigma_s(1, 1, 0, 5.65042e-02)
mat.set_sigma_s(1, 1, 1, 1.62452e-00)
# 2. Control Rod (Widthdrawn = Void)
mat.set_sigma_t(2, 0,    1.28407e-02)
mat.set_sigma_t(2, 1,    1.20676e-02)
mat.set_sigma_s(2, 0, 0, 1.27700e-02)
mat.set_sigma_s(2, 1, 0, 2.40997e-05)
mat.set_sigma_s(2, 1, 1, 1.07387e-02)
#
mat.finalize()

try:
  io = IO_HDF5("takeda_withdrawn.h5")
  io.open()
  io.write(inp)
  io.write(mat)
  io.write(mesh)
  io.close()
except:
  print("Error using HDF5---perhaps detran is not built with it."

# Solve
start = time.time()
solver = Eigen3D(inp, mat, mesh)
solver.solve()
print "elapsed = ", time.time() - start
state = solver.state()

try :
  silo = SiloOutput(mesh)
  silo.initialize("takeda_withdrawn.silo")
  silo.write_mesh_map("MATERIAL")
  silo.write_scalar_flux(state)
  silo.finalize()
except :
  print "Error using Silo---perhaps detran is not built with it."

#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
Manager.finalize()

