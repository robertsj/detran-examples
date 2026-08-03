# pyexamples/run_pin.py
#
# Solves a fully reflected C5G7 pin eigenproblem using diffusion.
#

import numpy as np
import time
from detran import *
try:
  from .pins_c5g7 import get_pins
  from .material_c5g7 import get_materials
  from .plot_utils import *
except ImportError:
  from pins_c5g7 import get_pins
  from material_c5g7 import get_materials
  from plot_utils import *

def run():
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
  problem = DiffusionEigenvalue(
    mesh=mesh,
    materials=mat,
    boundary="reflect",
    solver=EigenSolver.arnoldi(max_iters=1000, tolerance=1e-9),
  )
  result = problem.solve()
  print("elapsed = ", time.time() - start)

  #-----------------------------------------------------------------------------#
  # Plot
  #-----------------------------------------------------------------------------#
  try :
    state = result.state
    silo = SiloOutput(mesh)
    silo.initialize("c5g7pin.silo")
    silo.write_scalar_flux(state)
    silo.finalize()
  except :
    print("Silo error?")
  return result.keff


if __name__ == "__main__":
  run()
