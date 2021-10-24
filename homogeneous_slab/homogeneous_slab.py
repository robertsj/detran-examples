"""
Here, a homogeneous, one-group slab is defined and analyzed using
a few different techniques.
"""

from detran import *
import time

def get_input():
    """ fills default input parameters.
    """

    # Input
    inp = InputDB.Create()
    inp.put_str("problem_type",       "fixed")
    inp.put_int("number_groups",      1)
    inp.put_str("equation",           "dd")
    inp.put_str("inner_solver",       "SI")
    inp.put_int("inner_max_iters",    10000)
    inp.put_dbl("inner_tolerance",    1e-12)
    inp.put_int("inner_use_pc",       1)
    inp.put_int("inner_print_level",  0)
    # Note, we print out the outer, which gives the sweep count.  For
    # 1 group, though, there is no real outer iteration.
    inp.put_int("outer_print_out",            1)
    inp.put_str("quad_type",                  "gl")
    inp.put_int("quad_number_polar_octant",   1)
    inp.put_str("bc_west",                    "vacuum")
    #
    db = InputDB.Create("callow_db")
    db.put_dbl("linear_solver_atol",              1e-9);
    db.put_dbl("linear_solver_rtol",              1e-9);
    db.put_str("linear_solver_type",              "petsc");
    db.put_int("linear_solver_maxit",             5000);
    db.put_int("linear_solver_gmres_restart",     30);
    db.put_int("linear_solver_monitor_level",     0);
    db.put_str("pc_type",                         "petsc_pc");
    db.put_str("petsc_pc_type",                   "lu");
    inp.put_spdb("inner_solver_db",               db)
    inp.put_spdb("inner_pc_db",                   db)
    inp.put_spdb("outer_solver_db",               db)
    inp.put_spdb("eigen_solver_db",               db)
    return inp

def get_mesh():
    # Mesh
    mesh = Mesh1D.Create([200], [0.0, 100.0], [0])
    return mesh

def get_material(c=0.5):
    mat = Material.Create(1, 1)
    mat.set_sigma_t(0, 0, 1.0)
    mat.set_sigma_s(0, 0, 0, c)
    mat.finalize()
    return mat

def run():
    inp = get_input()
    mat = get_material()
    mesh = get_mesh()
    solver = Fixed1D(inp, mat, mesh)
    solver.setup()
    quad = solver.quadrature()
    q_e = ConstantSource.Create(1, mesh, 1.0, quad)
    solver.set_source(q_e)
    solver.set_solver()
    start = time.time()
    solver.solve()
    state = solver.state()
    elapsed = (time.time() - start)
    print elapsed, " seconds"

    plot_mesh_function(mesh, solver.state().phi(0))


if __name__ == "__main__":
  Manager.initialize(sys.argv)
  run()
