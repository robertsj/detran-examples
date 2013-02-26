# detran-examples/diffusion_benchmarks/lra2d.py
#
# Solves 2D LRA benchmark.
#
# Reference keff ~ 0.99636

from detran import *
import time

def run():

  #-----------------------------------------------------------------------------#
  # Input
  #-----------------------------------------------------------------------------#
  inp = InputDB.Create()
  inp.put_int("number_groups",                      2)
  inp.put_int("dimension",                          2)
  inp.put_str("equation",                           "diffusion")
  inp.put_str("bc_west",                            "reflect")
  inp.put_str("bc_east",                            "vacuum")
  inp.put_str("bc_south",                           "reflect")
  inp.put_str("bc_north",                           "vacuum")
  inp.put_int("eigen_max_iters",                    1000)
  inp.put_str("eigen_solver",                       "arnoldi")
  db = InputDB.Create("callow_db")
  # outer gmres parameters
  db.put_dbl("linear_solver_atol",                  1e-8);
  db.put_dbl("linear_solver_rtol",                  1e-8);
  db.put_str("linear_solver_type",                  "petsc");
  db.put_int("linear_solver_maxit",                 5000);
  db.put_int("linear_solver_gmres_restart",         30);
  db.put_str("eigen_solver_type",                   "slepc");
  db.put_int("linear_solver_monitor_level",         0);
  inp.put_spdb("outer_solver_db", db)
  inp.put_spdb("eigen_solver_db", db)

  #-----------------------------------------------------------------------------#
  # Material
  #-----------------------------------------------------------------------------#
  # Note, because this problem has the removal cross section specified, but
  # not any within-group scattering, we'll simply put the removal values into
  # the total, since removal is total minus within group.  Note also that the
  # removal values are adjusted for axial buckling.
  mat = Material.Create(5, 2, "LRA-2D")
  B = 1.e-4
  # Material 0
  mat.set_sigma_t(0,   vec_dbl([0.008252, 0.1003]))
  mat.set_sigma_s(0, 1, 0,      0.02533) 
  mat.set_diff_coef(0, vec_dbl([1.255, 0.211]))
  mat.set_sigma_f(0,   vec_dbl([0.004602, 0.1091 ]))
  mat.set_chi(0, 0, 1.0)
  # Material 1
  mat.set_sigma_t(1,   vec_dbl([0.007181, 0.07047]))
  mat.set_sigma_s(1, 1, 0,      0.02767) 
  mat.set_diff_coef(1, vec_dbl([1.268, 0.1902]))
  mat.set_sigma_f(1,   vec_dbl([0.004609, 0.08675]))
  mat.set_chi(1, 0, 1.0)
  # Material 2
  mat.set_sigma_t(2,   vec_dbl([0.008002, 0.08344]))
  mat.set_sigma_s(2, 1, 0,      0.02617) 
  mat.set_diff_coef(2, vec_dbl([1.259, 0.2091]))
  mat.set_sigma_f(2,   vec_dbl([0.004663, 0.1021]))
  mat.set_chi(2, 0, 1.0)
  # Material 3
  mat.set_sigma_t(3,   vec_dbl([0.008002, 0.073324 ]))
  mat.set_sigma_s(3, 1, 0,      0.02617) 
  mat.set_diff_coef(3, vec_dbl([1.259, 0.2091]))
  mat.set_sigma_f(3,   vec_dbl([0.004663, 0.1021 ]))
  mat.set_chi(3, 0,             1.0)
  # Material 4
  mat.set_sigma_t(4,   vec_dbl([0.0006034, 0.01911]))
  mat.set_sigma_s(4, 1, 0,      0.04754) 
  mat.set_diff_coef(4, vec_dbl([1.257, 0.1592]))
  # Corrections to removal
  for m in range(0, 5) :
    # group 0 -- buckling
    orig = mat.sigma_t(m, 0)
    new  = orig + mat.diff_coef(m, 0)*B
    mat.set_sigma_t(m, 0, new)
    # group 0 -- s12 
    orig = mat.sigma_t(m, 0)
    new  = orig + mat.sigma_s(m, 1, 0)
    mat.set_sigma_t(m, 0, new)
    # group 1 -- buckling
    orig = mat.sigma_t(m, 1)
    new  = orig + mat.diff_coef(m, 1)*B
    mat.set_sigma_t(m, 1, new)
  mat.finalize()

  #-----------------------------------------------------------------------------#
  # Geometry
  #-----------------------------------------------------------------------------#
  cm = [0.0, 15.0, 75.0, 105.0, 120.0, 135.0, 165.0]
  fm = vec_int(6, 10)
  for i in range(1, 6) :
    fm[i] = 2*fm[i]
  mt = [1, 0, 1, 2, 2, 4,
        0, 0, 0, 2, 2, 4,
        1, 0, 1, 2, 2, 4,
        2, 2, 2, 3, 4, 4,
        2, 2, 2, 4, 4, 4,
        4, 4, 4, 4, 4, 4]
  mesh = Mesh2D.Create(fm, fm, cm, cm, mt)

  #-----------------------------------------------------------------------------#
  # Execute
  #-----------------------------------------------------------------------------#
  Manager.initialize(sys.argv)
  solver = Eigen2D(inp, mat, mesh)
  t = time.time()
  solver.solve()
  state = solver.state()
  print "elapsed = ", time.time()-t

  #-----------------------------------------------------------------------------#
  # Plot
  #-----------------------------------------------------------------------------#
  try :
    silo = SiloOutput(mesh)
    silo.initialize("lra2d.silo")
    silo.write_scalar_flux(state)
    silo.finalize()
  except :
    print "Silo error (not installed?)"

if __name__ == "__main__":
  Manager.initialize(sys.argv)
  run()
