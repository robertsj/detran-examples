# detran-examples/diffusion_benchmarks/iaea.py
#
# Solves 2D IAEA benchmark.
#
# Reference keff ~ 1.02959

import numpy as np
import time
import sys
from detran import *

#-----------------------------------------------------------------------------#
# Input
#-----------------------------------------------------------------------------#
inp = InputDB.Create()
inp.put_int("number_groups",            2)
inp.put_int("diffusion_max_iters",      1000)
inp.put_dbl("diffusion_tolerance",      1e-6)
inp.put_str("bc_west",                  "reflect")
inp.put_str("bc_east",                  "vacuum")
inp.put_str("bc_south",                 "reflect")
inp.put_str("bc_north",                 "vacuum")

#-----------------------------------------------------------------------------#
# Material
#-----------------------------------------------------------------------------#
# Note, because this problem has the removal cross section specified, but
# not any within-group scattering, we'll simply put the removal values into
# the total, since removal is total minus within group.  Note also that the
# removal values are adjusted for axial buckling.
mat = Material.Create(4, 2, False)
B = 1.e-4
# Material 0
mat.set_sigma_t(0,   vec_dbl([0.008252, 0.1003]))
mat.set_sigma_s(0, 1, 0,      0.02533) 
mat.set_diff_coef(0, vec_dbl([1.255, 0.211]))
mat.set_sigma_f(0,   vec_dbl([0.004602, 0.1091 ]))
mat.set_chi(0, 0, 1.0)
# Material 1
mat.set_sigma_t(1,   vec_dbl([0.007181, 0.07047]))
mat.set_sigma_s(1, 1, 0,      0.02767) 
mat.set_diff_coef(1, vec_dbl([1.268, 0.1902]))
mat.set_sigma_f(1,   vec_dbl([0.004609, 0.08675]))
mat.set_chi(1, 0, 1.0)
# Material 2
mat.set_sigma_t(2,   vec_dbl([0.008002, 0.08344]))
mat.set_sigma_s(2, 1, 0,      0.02617) 
mat.set_diff_coef(2, vec_dbl([1.259, 0.2091]))
mat.set_sigma_f(2,   vec_dbl([0.004663, 0.1021]))
mat.set_chi(2, 0, 1.0)
# Material 3
mat.set_sigma_t(3,   vec_dbl([0.008002, 0.073324 ]))
mat.set_sigma_s(3, 1, 0,      0.02617) 
mat.set_diff_coef(3, vec_dbl([1.259, 0.2091]))
mat.set_sigma_f(3,   vec_dbl([0.004663, 0.1021 ]))
mat.set_chi(3, 0,             1.0)
# Material 4
mat.set_sigma_t(4,   vec_dbl([0.0006034, 0.01911]))
mat.set_sigma_s(4, 1, 0,      0.04754) 
mat.set_diff_coef(4, vec_dbl([1.257, 0.1592]))
# correct removal for buckling and 1->2 scatter
for m in range(0, 5) :
  orig = mat.sigma_t(m, 1)
  new  = orig + mat.diff_coef(m, 0)*B
  mat.set_sigma_t(m, 1, new)
  orig = mat.sigma_t(m, 0)
  new  = orig + mat.sigma_s(m, 1, 0)
  mat.set_sigma_t(m, 0, new)
mat.finalize()

#-----------------------------------------------------------------------------#
# Geometry
#-----------------------------------------------------------------------------#
cm = [0.0, 15.0, 75.0, 105.0, 120.0, 135.0, 165.0]
fm = vec_int(6, 20)
for i in range(1, 6) :
  fm[i] = 2*fm[i]
mt = [1, 0, 1, 2, 2, 4,
      0, 0, 0, 2, 2, 4,
      1, 0, 1, 2, 2, 4,
      2, 2, 2, 3, 4, 4,
      2, 2, 2, 4, 4, 4,
      4, 4, 4, 4, 4, 4]
mesh = Mesh2D.Create(fm, fm, cm, cm, mt)

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
plot_mesh_function(mesh, np.asarray(state.phi(1)))
exit()
silo = SiloOutput(mesh)
silo.initialize("lra2d.silo")
silo.write_scalar_flux(state)
silo.finalize()

#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
Manager.finalize()

