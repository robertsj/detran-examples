# detran-examples/diffusion_benchmarks/iaea2d.py
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
inp.put_int("number_groups",              2)
inp.put_int("diffusion_max_iters",        1000)
inp.put_dbl("diffusion_tolerance",        1e-6)
inp.put_str("diffusion_eigensolver_type", "slepc")
inp.put_str("bc_west",                    "reflect")
inp.put_str("bc_east",                    "vacuum")
inp.put_str("bc_south",                   "reflect")
inp.put_str("bc_north",                   "vacuum")

#-----------------------------------------------------------------------------#
# Material
#-----------------------------------------------------------------------------#
# Note, because this problem has the removal cross section specified, but
# not any within-group scattering, we'll simply put the removal values into
# the total, since removal is total minus within group.  Note also that the
# removal values are adjusted for axial buckling.
mat = Material.Create(4, 2, False)
# Material 0
mat.set_sigma_t(0, vec_dbl([0.030120, 0.080032]))
mat.set_sigma_s(0, 1, 0, 0.02) 
mat.set_diff_coef(0, vec_dbl([1.5, 0.4]))
mat.set_sigma_f(0, 1, 0.135)
mat.set_chi(0, 0, 1.0)
# Material 1
mat.set_sigma_t(1, vec_dbl([0.030120, 0.085032]))
mat.set_sigma_s(1, 1, 0, 0.02) 
mat.set_diff_coef(1, vec_dbl([1.5, 0.4]))
mat.set_sigma_f(1, 1, 0.135)
mat.set_chi(1, 0, 1.0)
# Material 2
mat.set_sigma_t(2, vec_dbl([0.030120, 0.130032]))
mat.set_sigma_s(2, 1, 0, 0.02) 
mat.set_diff_coef(2, vec_dbl([1.5, 0.4]))
mat.set_sigma_f(2, 1, 0.135)
mat.set_chi(2, 0, 1.0)
## Material 3
mat.set_sigma_t(3, vec_dbl([0.040160, 0.010024]))
mat.set_sigma_s(3, 1, 0, 0.04) 
mat.set_diff_coef(3, vec_dbl([2.0, 0.3]))
mat.finalize()

#-----------------------------------------------------------------------------#
# Geometry
#-----------------------------------------------------------------------------#
cm = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0]
fm = vec_int(9, 5)
for i in range(1, 9) :
  fm[i] = 2*fm[i]
mt = [2, 1, 1, 1, 2, 1, 1, 0, 3,
      1, 1, 1, 1, 1, 1, 1, 0, 3,
      1, 1, 1, 1, 1, 1, 0, 0, 3,
      1, 1, 1, 1, 1, 1, 0, 3, 3,
      2, 1, 1, 1, 2, 0, 0, 3, 3,
      1, 1, 1, 1, 0, 0, 3, 3, 3,
      1, 1, 0, 0, 0, 3, 3, 3, 3,
      0, 0, 0, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3, 3]
mesh = Mesh2D.Create(fm, fm, cm, cm, mt)

#-----------------------------------------------------------------------------#
# State
#-----------------------------------------------------------------------------#
state = State.Create(inp, mesh)

#-----------------------------------------------------------------------------#
# Execute
#-----------------------------------------------------------------------------#
Manager.initialize(sys.argv)
solver = DiffusionEigen2D(inp, mat, mesh, state)
t = time.time()
solver.solve()
print "elapsed = ", time.time()-t

#-----------------------------------------------------------------------------#
# Plot
#-----------------------------------------------------------------------------#
plot_mesh_function(mesh, np.asarray(state.phi(1)))
silo = SiloOutput(mesh)
silo.initialize("iaea2d.silo")
silo.write_scalar_flux(state)
silo.finalize()

#-----------------------------------------------------------------------------#
# Wrap Up
#-----------------------------------------------------------------------------#
Manager.finalize()

