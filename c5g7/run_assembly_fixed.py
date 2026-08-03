# pyexamples/run_assembly.py
#
# Test the efficacy of preconditioners for multigroup fixed 
# source multiplying problems.

import numpy as np
import time
from detran import *
try:
  from .assemblies_c5g7 import get_assemblies
  from .pins_c5g7 import get_pins
  from .material_c5g7 import get_materials
except ImportError:
  from assemblies_c5g7 import get_assemblies
  from pins_c5g7 import get_pins
  from material_c5g7 import get_materials
import cProfile

def run() :

  #for n in range(1, 11) :
    n = 2
    #-----------------------------------------------------------------------------#
    # Input 
    #-----------------------------------------------------------------------------#
    solver_db = SolverDatabase(
      linear=LinearSolver.gmres(
        atol=0.0,
        rtol=1e-8,
        max_iters=50000,
        restart=30,
        monitor_level=1,
      )
    )
    preconditioner_db = SolverDatabase(
      linear=LinearSolver.gmres(
        atol=0.0,
        rtol=0.1,
        max_iters=5000,
        restart=30,
        monitor_level=0,
        pc_type="ilut",
        pc_side=1,
        ilu_levels=2,
        ilu_thresh=0.0,
      )
    )
    inp = fixed_source_input(
      groups=7,
      dimension=2,
      equation="dd",
      boundary="vacuum",
      inner=IterationSolver.source_iteration(
        max_iters=1000000,
        tolerance=1e-8,
        print_level=1,
        print_interval=10,
        pc_type="DSA",
        pc_side=2,
      ),
      outer=IterationSolver.gauss_seidel(
        max_iters=1000,
        tolerance=1e-8,
        print_level=1,
        print_interval=1,
      ),
      quadrature=Quadrature(polar=n, azimuth=n),
      solver_databases={
        "inner_solver_db": solver_db,
        "inner_pc_db": preconditioner_db,
        "outer_solver_db": solver_db,
        "outer_pc_db": preconditioner_db,
      },
    )
    #-----------------------------------------------------------------------------#
    # Material
    #-----------------------------------------------------------------------------#
    mat = get_materials()
    mat.compute_diff_coef()
    #-----------------------------------------------------------------------------#
    # Geometry
    #-----------------------------------------------------------------------------#
    assemblies = get_assemblies(7, True)
    mesh = assemblies[0].mesh()
    #-----------------------------------------------------------------------------#
    # Source and Solve
    #-----------------------------------------------------------------------------#
    solver = Fixed2D(inp, mat, mesh)
    solver.setup()
    quad = solver.quadrature()
    q_e = ConstantSource(7, mesh, 1.0, quad)
    solver.set_source(q_e)
    solver.set_solver()
    t = time.time()
    solver.solve()
    print("elapsed = ", time.time()-t)
    phi = np.asarray(solver.state().phi(0))
    #plot_mesh_function(mesh, phi)

if __name__ == "__main__":
  #Manager.initialize(sys.argv)
  run()
