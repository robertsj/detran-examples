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
from detran import *
try:
  from .assemblies_c5g7 import get_assemblies
  from .pins_c5g7 import get_pins
  from .material_c5g7 import get_materials
  from .plot_utils import *
except ImportError:
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
  assemblies = get_assemblies(7, True)
  mesh = assemblies[0].mesh()
  #-----------------------------------------------------------------------------#
  # Execute
  #-----------------------------------------------------------------------------#
  #Manager.initialize(sys.argv)
  problem = DiffusionEigenvalue(
    mesh=mesh,
    materials=mat,
    boundary="reflect",
    solver=EigenSolver.arnoldi(
      max_iters=1000,
      inner=LinearSolver.gmres(
        atol=1e-8,
        rtol=1e-8,
        max_iters=5000,
        restart=30,
      ),
    ),
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
  return result.keff


if __name__ == "__main__":
  run()
