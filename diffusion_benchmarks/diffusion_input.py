"""
Base input database for diffusion problems.
"""

from detran import InputDB

def get_input(solver=None):
    """ Default parameters for solving quarter-core, diffusion eigenproblems.
    """

    inp = InputDB.Create()
    inp.put_int("number_groups",                      2)
    inp.put_int("dimension",                          2)
    inp.put_str("equation",                           "diffusion")
    inp.put_str("bc_west",                            "reflect")
    inp.put_str("bc_east",                            "vacuum")
    inp.put_str("bc_south",                           "reflect")
    inp.put_str("bc_north",                           "vacuum")
    inp.put_int("eigen_max_iters",                    1000)
    ## TODO: callow
    inp.put_str("eigen_solver",                       "arnoldi")
    db = InputDB.Create("callow_db")
    # outer gmres parameters
    db.put_dbl("linear_solver_atol",                  1e-8);
    db.put_dbl("linear_solver_rtol",                  1e-8);
    db.put_str("linear_solver_type",                  "gmres");
    db.put_int("linear_solver_maxit",                 5000);
    db.put_int("linear_solver_gmres_restart",         30);
    db.put_str("eigen_solver_type",                   "gd");
    db.put_int("eigen_solver_monitor_level",          2);
    db.put_int("linear_solver_monitor_level",         0);
    inp.put_spdb("outer_solver_db", db)
    inp.put_spdb("eigen_solver_db", db)

    if solver == "slepc":
        pass

    return inp




if __name__ == "__main__":
    inp = get_input()
    inp.display()
