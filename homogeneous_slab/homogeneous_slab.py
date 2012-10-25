# pyexamples/homogeneous_slab/homogeneous_slab.py
#
# This problem defines a homogeneous slab with varying scattering 
# ratios in the one group approximation.  It allows us to test 
# various fixed source solvers (e.g. SI and GMRES) along with
# preconditioning.

import numpy as np
import time
import sys
from detran import *

# Initialize
Manager.initialize(sys.argv)

# Input
inp = InputDB.Create()
inp.put_str("problem_type",       "fixed")
inp.put_int("number_groups",      1)
inp.put_str("equation",           "dd")
inp.put_str("inner_solver",       "SI")
inp.put_int("inner_max_iters",    10000)
inp.put_dbl("inner_tolerance",    1e-8)
inp.put_int("inner_use_pc",       1)
inp.put_int("inner_print_level",  0)
# Note, we print out the outer, which gives the sweep count.  For 
# 1 group, though, there is no real outer iteration.
inp.put_int("outer_print_out",            1)
inp.put_str("quad_type",                  "gausslegendre")
inp.put_int("quad_number_polar_octant",   1)
inp.put_str("bc_west",                    "reflect")

# Mesh
mesh = Mesh1D.Create([200], [0.0, 100.0], [0])

# Loop over scattering ratios
c = np.linspace(0.0, 0.99, 10)
for i in range(0, 10) :

  # Material
  mat = Material.Create(1, 1, False)
  mat.set_sigma_t(0, 0, 1.0)
  mat.set_sigma_s(0, 0, 0, c[i])
  mat.finalize()

  # Initialize and set source
  solver = Fixed1D(inp, mat, mesh)
  solver.setup()
  quad = solver.quadrature()
  q_e = ConstantSource.Create(1, mesh, 1.0, quad)
  solver.set_source(q_e)
  solver.set_solver()

  # Solve
  start = time.time()
  solver.solve()
  elapsed = (time.time() - start)
  print elapsed, " seconds"

Manager.finalize()


