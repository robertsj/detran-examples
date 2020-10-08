# pyexamples/run_core.py
#
# Solves the 2-D C5G7 benchmark problem.
#
# Using 7x7 volume-conserving pin cell, DD SN, and a 3x6 QR quadrature,
# the reference keff is 1.186572165.

import numpy as np
import time
import sys
from detran import *
from core_c5g7 import get_core
from assemblies_c5g7 import get_assemblies
from pins_c5g7 import get_pins
from material_c5g7 import get_materials

def run() :

    inp = InputDB.Create()
    inp.put_str("equation",                       "dd")
    inp.put_str("problem_type",                   "eigenvalue")
    inp.put_int("number_groups",                  7)
    #
    inp.put_str("inner_solver",                   "GMRES")
    inp.put_int("inner_max_iters",                1)
    inp.put_dbl("inner_tolerance",                1e-3)
    inp.put_int("inner_print_level",              0)
    inp.put_int("inner_print_interval",           10)
    #
    inp.put_str("outer_solver",                   "GMRES")
    inp.put_int("outer_max_iters",                1000)
    inp.put_dbl("outer_tolerance",                1e-8)
    inp.put_int("outer_print_level",              0)
    inp.put_int("outer_print_interval",           1)
    inp.put_int("outer_krylov_group_cutoff",      0) # gmres on whole mg
    inp.put_str("outer_pc_type",                  "mgcmdsa")

    inp.put_str("eigen_solver_pc_type",                  "mgcmdsa")
    inp.put_int('mgpc_coarse_mesh_level',         7)
    inp.put_int('mgpc_condensation_option',       1)
    inp.put_int('mgpc_cmdsa_use_smoothing',       1)
    inp.put_int('mgpc_cmdsa_smoothing_iters',     3)
    inp.put_dbl('mgpc_cmdsa_smoothing_relax',     0.7)

    #
    inp.put_str("eigen_solver",                   "PI")
    inp.put_int("eigen_max_iters",                1000)
    inp.put_dbl("eigen_tolerance",                1e-6)
    inp.put_int("eigen_print_level",              2)
    inp.put_int("eigen_print_interval",           2)
    #
    inp.put_str("bc_west",                        "reflect")
    inp.put_str("bc_east",                        "vacuum")
    inp.put_str("bc_south",                       "reflect")
    inp.put_str("bc_north",                       "vacuum")
    #
    inp.put_int("quad_number_polar_octant",       3)
    inp.put_int("quad_number_azimuth_octant",     6)
    #
    db = InputDB.Create("callow_db")
    db.put_dbl("linear_solver_atol",              1e-9);
    db.put_dbl("linear_solver_rtol",              1e-9);
    db.put_str("linear_solver_type",              "gmres");
    db.put_int("linear_solver_maxit",             5000);
    db.put_int("linear_solver_gmres_restart",     30);
    db.put_int("linear_solver_monitor_level",     1);
    db.put_str("pc_type",                         "ilu0");
    #db.put_str("petsc_pc_type",                   "lu");
    db.put_str("eigen_solver_type",               "gd");
    db.put_int("eigen_solver_monitor_level",      1);
    inp.put_spdb("inner_solver_db",               db)
    inp.put_spdb("inner_pc_db",                   db)
    inp.put_spdb("outer_solver_db",               db)
    inp.put_spdb("outer_pc_db",               db)
    inp.put_spdb("eigen_solver_db",               db)
    inp.put_spdb("eigen_solver_pc_db",            db)



    mat = get_materials()

    core = get_core(7, True)
    mesh = core.mesh()

    start = time.time()
    solver = Eigen2D(inp, mat, mesh)
    solver.solve()
    print "elapsed = ", time.time() - start


if __name__ == "__main__":
  Manager.initialize(sys.argv)
  run()
  Manager.finalize()
