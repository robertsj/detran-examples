"""
Consolidate input generation.
"""

from detran import *
import pickle

from slab_materials import *

def get_input(case='diffusion'):
    inp = get_base_input()
    if case == 'diffusion':
        inp.put_str("equation", "diffusion")
    elif case == 'transport_full':
        inp.put_str("equation", "dd")
        inp.put_int("ts_discrete", 1)
        inp.put_int("store_angular_flux", 1)

    elif case == 'transport_iso':
        inp.put_str("equation", "dd")
        inp.put_int("ts_discrete", 0)
        inp.put_int("store_angular_flux", 0)
    return inp

def get_base_input() :
    inp = utilities.InputDB.Create()
    inp.put_int("number_groups",                  2)
    inp.put_int("dimension",                      1)
    inp.put_str("equation",                       "diffusion")
    inp.put_int("adjoint",                        0)
    inp.put_str("bc_west",                        "vacuum")
    inp.put_str("bc_east",                        "vacuum")
    inp.put_int("quad_number_polar_octant",       8)
    inp.put_str("eigen_solver",                   "PI")
    inp.put_dbl("eigen_tolerance",                1e-14)
    inp.put_int("eigen_max_iters",                1000)
    inp.put_str("outer_solver",                   "GMRES")
    inp.put_dbl("outer_tolerance",                1e-12)
    inp.put_int("outer_max_iters",                1000)
    inp.put_int("outer_print_level",              0)
    inp.put_int("outer_krylov_group_cutoff",      0)
    inp.put_str("outer_pc_type",                  "none")
    inp.put_int('mgpc_coarse_mesh_level',         5)
    inp.put_int('mgpc_condensation_option',       1)
    inp.put_int('mgpc_cmdsa_use_smoothing',       1)
    inp.put_int('mgpc_cmdsa_smoothing_iters',     2)
    inp.put_dbl('mgpc_cmdsa_smoothing_relax',     0.8)
    inp.put_str("inner_solver",                   "GMRES")
    inp.put_dbl("inner_tolerance",                1e-12)
    inp.put_int("inner_max_iters",                1000)
    inp.put_int("inner_print_level",              0)
    # gmres parameters
    db = utilities.InputDB.Create("callow_db")
    db.put_dbl("linear_solver_atol",              0.0)
    db.put_dbl("linear_solver_rtol",              1e-13)
    db.put_str("linear_solver_type",              "gmres")
    db.put_int("linear_solver_maxit",             1000)
    db.put_int("linear_solver_gmres_restart",     50)
    db.put_int("linear_solver_monitor_level",     0)
    db.put_str("eigen_solver_type",               "power")
    db.put_int("eigen_solver_maxit",              1000)
    db.put_int("eigen_solver_monitor_level",      2)
    db.put_dbl("eigen_solver_tol",                1.0e-14)
    inp.put_spdb("inner_solver_db", db)
    inp.put_spdb("inner_pc_db", db)
    inp.put_spdb("outer_solver_db", db) 
    inp.put_spdb("eigen_solver_db", db)
    inp.put_int("ts_max_steps",                   1000000)
    inp.put_int("ts_scheme",                      Time1D.BDF1)
    inp.put_int("ts_output",                      0)
    inp.put_dbl("ts_step_size",                   0.05)
    inp.put_dbl("ts_final_time",                  60.0)
    inp.put_int("ts_max_iters",                   10)
    inp.put_dbl("ts_tolerance",                   1.0e-5)
    inp.put_int("ts_discrete",                    0)
    inp.put_int("store_angular_flux",             0)
    
    return inp