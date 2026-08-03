# pyexamples/run_core.py
#
# Solves the 2-D C5G7 benchmark problem via diffusion.
#
# Using 7x7 volume-conserving pin cell, DD SN, and a QR quadrature
# of order 18, the reference keff is 1.186572165.

import numpy as np
import time
from detran import *
try:
  from .core_c5g7 import get_core
  from .assemblies_c5g7 import get_assemblies
  from .pins_c5g7 import get_pins
  from .material_c5g7 import get_materials
  from .plot_utils import *
except ImportError:
  from core_c5g7 import get_core
  from assemblies_c5g7 import get_assemblies
  from pins_c5g7 import get_pins
  from material_c5g7 import get_materials
  from plot_utils import *
def run():
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
  core = get_core(7, True)
  mesh = core.mesh()
  #-----------------------------------------------------------------------------#
  # Execute
  #-----------------------------------------------------------------------------#
  #Manager.initialize(sys.argv)
  problem = DiffusionEigenvalue(
    mesh=mesh,
    materials=mat,
    boundary=("reflect", "vacuum", "reflect", "vacuum"),
    solver=EigenSolver.arnoldi(max_iters=1000, tolerance=1e-6),
  )
  t = time.time()
  result = problem.solve()
  print("elapsed = ", time.time()-t)
  #-----------------------------------------------------------------------------#
  # Plot
  #-----------------------------------------------------------------------------#
  state = result.state
  plot_flux(mesh, np.asarray(state.phi(0)))
  #-----------------------------------------------------------------------------#
  # Wrap Up
  #-----------------------------------------------------------------------------#
  #Manager.finalize()
  return state


if __name__ == "__main__":
  run()
