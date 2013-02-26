# pyexamples/run_assembly.py
#
# Test the efficacy of preconditioners for multigroup fixed 
# source multiplying problems.

import numpy as np
import time
import sys
from detran import *
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
    inp = utilities.InputDB.Create()
    inp.put_str("equation",                   "dd")
    inp.put_str("problem_type",               "eigenvalue")
    inp.put_int("number_groups",              7)
    #
    inp.put_str("inner_solver",               "SI")
    inp.put_int("inner_max_iters",            1000000)
    inp.put_dbl("inner_tolerance",            1e-8)
    inp.put_int("inner_print_level",          1)
    inp.put_int("inner_print_interval",       10)
    inp.put_str("inner_pc_type",              "DSA")
    inp.put_int("inner_pc_side",              2)
    #
    inp.put_str("outer_solver",               "GS")
    inp.put_int("outer_max_iters",            1000)
    inp.put_dbl("outer_tolerance",            1e-8)
    inp.put_int("outer_print_level",          1)
    inp.put_int("outer_print_interval",       1)
    #inp.put_str("outer_pc_type",              "mgdsa")
    #inp.put_int("outer_krylov_group_cutoff",  0)
    #inp.put_int("outer_pc_side",              2)
    #
    inp.put_str("bc_west",                    "vacuum")
    inp.put_str("bc_east",                    "vacuum")
    inp.put_str("bc_south",                   "vacuum")
    inp.put_str("bc_north",                   "vacuum")
    #
    inp.put_str("quad_type",                  "chebyshevlegendre")
    inp.put_int("quad_number_polar_octant",   n)
    inp.put_int("quad_number_azimuth_octant", n)
    #
    solver_db = utilities.InputDB.Create("solver_db")
    solver_db.put_dbl("linear_solver_atol",              0.0);
    solver_db.put_dbl("linear_solver_rtol",              1e-8);
    solver_db.put_str("linear_solver_type",              "petsc");
    solver_db.put_int("linear_solver_maxit",             50000);
    solver_db.put_int("linear_solver_gmres_restart",     30);
    solver_db.put_int("linear_solver_monitor_level",     1);
    #
    preconditioner_db = utilities.InputDB.Create("preconditioner_db")
    preconditioner_db.put_dbl("linear_solver_atol",              0.0);
    preconditioner_db.put_dbl("linear_solver_rtol",              0.1);
    preconditioner_db.put_str("linear_solver_type",              "petsc");
    preconditioner_db.put_int("linear_solver_maxit",             5000);
    preconditioner_db.put_int("linear_solver_gmres_restart",     30);
    preconditioner_db.put_int("linear_solver_monitor_level",     0);
    preconditioner_db.put_str("pc_type",                         "petsc_pc");
    preconditioner_db.put_str("petsc_pc_type",                   "ilu");
    preconditioner_db.put_int("petsc_pc_factor_levels",          2);
    #
    inp.put_spdb("inner_solver_db", solver_db)
    inp.put_spdb("inner_pc_db",     preconditioner_db)
    inp.put_spdb("outer_solver_db", solver_db)
    inp.put_spdb("outer_pc_db",     preconditioner_db)
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
    q_e = external_source.ConstantSource.Create(7, mesh, 1.0, quad)
    solver.set_source(q_e)
    solver.set_solver()
    t = time.time()
    solver.solve()
    print "elapsed = ", time.time()-t
    phi = np.asarray(solver.state().phi(0))
    #plot_mesh_function(mesh, phi)

if __name__ == "__main__":
  Manager.initialize(sys.argv)
  run()
