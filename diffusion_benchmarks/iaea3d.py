# detran-examplesdiffusion_benchmarks/iaea.py
#
# Solves the 3D IAEA PWR benchmark problem.  Note, because of modeling
# limitations, the void regions beyond the normal reflector are simply
# filled with more reflector.
#
# Reference keff ~ 1.029096
#
# Reference: 
#   Benchmark Problem Book, ANL-7416, Suppl. 2, Argonne National 
#   Laboratory (1977)

from detran import *
import time

def run():

  #-----------------------------------------------------------------------------#
  # Input
  #-----------------------------------------------------------------------------#

  inp = InputDB.Create()
  inp.put_int("number_groups",                  2)
  inp.put_int("dimension",                      2)
  inp.put_str("equation",                       "diffusion")
  inp.put_str("bc_west",                        "reflect")
  inp.put_str("bc_east",                        "reflect")
  inp.put_str("bc_south",                       "reflect")
  inp.put_str("bc_north",                       "reflect")
  inp.put_str("bc_bottom",                      "reflect")
  inp.put_str("bc_top",                         "reflect")
  inp.put_int("eigen_max_iters",                1000)
  inp.put_str("eigen_solver",                   "arnoldi")
  db = InputDB.Create("callow_db")
  # outer gmres parameters
  db.put_dbl("linear_solver_atol",                  1e-14);
  db.put_dbl("linear_solver_rtol",                  1e-14);
  db.put_str("linear_solver_type",                  "petsc");
  db.put_int("linear_solver_maxit",                 5000);
  db.put_int("linear_solver_gmres_restart",         30);
  db.put_str("eigen_solver_type",                   "slepc");
  db.put_dbl("eigen_solver_tol",                    1e-14);
  db.put_int("linear_solver_monitor_level",         0);
  inp.put_spdb("outer_solver_db", db)
  inp.put_spdb("eigen_solver_db", db)

  #-----------------------------------------------------------------------------#
  # Material
  #-----------------------------------------------------------------------------#

  # Note, all absorption cross sections are simply put into the total. 

  mat = Material.Create(6, 2, "IAEA-3D")
  # Fuel 1
  mat.set_sigma_t(0, vec_dbl([0.03, 0.08]))
  mat.set_sigma_s(0, 1, 0, 0.02) 
  mat.set_diff_coef(0, vec_dbl([1.5, 0.4]))
  mat.set_sigma_f(0, 1, 0.135)
  mat.set_chi(0, 0, 1.0)
  # Fuel 1 + Rod
  mat.set_sigma_t(1, vec_dbl([0.030, 0.085]))
  mat.set_sigma_s(1, 1, 0, 0.02) 
  mat.set_diff_coef(1, vec_dbl([1.5, 0.4]))
  mat.set_sigma_f(1, 1, 0.135)
  mat.set_chi(1, 0, 1.0)
  # Fuel 2
  mat.set_sigma_t(2, vec_dbl([0.03, 0.13]))
  mat.set_sigma_s(2, 1, 0, 0.02) 
  mat.set_diff_coef(2, vec_dbl([1.5, 0.4]))
  mat.set_sigma_f(2, 1, 0.135)
  mat.set_chi(2, 0, 1.0)
  # Reflector
  mat.set_sigma_t(3, vec_dbl([0.04, 0.01]))
  mat.set_sigma_s(3, 1, 0, 0.04) 
  mat.set_diff_coef(3, vec_dbl([2.0, 0.3]))
  # Reflector + Rod
  mat.set_sigma_t(4, vec_dbl([0.04, 0.55]))
  mat.set_sigma_s(4, 1, 0, 0.04) 
  mat.set_diff_coef(4, vec_dbl([2.0, 0.3]))
  # High Absorber
  mat.set_sigma_t(5, vec_dbl([1.00, 1.00]))
  mat.set_sigma_s(5, 1, 0, 0.00) 
  mat.set_diff_coef(5, vec_dbl([0.3333, 0.3333]))
  mat.finalize()

  #-----------------------------------------------------------------------------#
  # Geometry
  #-----------------------------------------------------------------------------#

  # This sets up for a 2cm mesh in all directions.

  # XY plane discretization
  cmH = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0]
  fmH = vec_int(9, 2)
  for i in range(1, 9) :
    fmH[i] = 2*fmH[i]
  # Axial discretization
  cmV = [0.0, 20.0, 280.0, 360.0, 380.0]
  #fmV = [  5,   65,    20,     5       ]
  fmV = [  2,   26,     8,     2       ]
  mt = [# 0.0 - 20.0
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        # 20.0 - 280.0
        2, 1, 1, 1, 2, 1, 1, 0, 3,
        1, 1, 1, 1, 1, 1, 1, 0, 3,
        1, 1, 1, 1, 1, 1, 0, 0, 3,
        1, 1, 1, 1, 1, 1, 0, 3, 3,
        2, 1, 1, 1, 2, 0, 0, 3, 3,
        1, 1, 1, 1, 0, 0, 3, 3, 3,
        1, 1, 0, 0, 0, 3, 3, 3, 3,
        0, 0, 0, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        # 280.0 - 360.0
        2, 1, 1, 1, 2, 1, 1, 0, 3,
        1, 1, 1, 1, 1, 1, 1, 0, 3,
        1, 1, 2, 1, 1, 1, 0, 0, 3,
        1, 1, 1, 1, 1, 1, 0, 3, 3,
        2, 1, 1, 1, 2, 0, 0, 3, 3,
        1, 1, 1, 1, 0, 0, 3, 3, 3,
        1, 1, 0, 0, 0, 3, 3, 3, 3,
        0, 0, 0, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        # 360.0 - 380.0
        4, 3, 3, 3, 4, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 4, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        4, 3, 3, 3, 4, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 3, 3, 3
       ]
  mesh = Mesh3D.Create(fmH, fmH, fmV, cmH, cmH, cmV, mt)

  #-----------------------------------------------------------------------------#
  # Execute
  #-----------------------------------------------------------------------------#
  solver = Eigen2D(inp, mat, mesh)
  t = time.time()
  solver.solve()
  state = solver.state()
  print state.eigenvalue()
  print "elapsed = ", time.time()-t

  #-----------------------------------------------------------------------------#
  # Plot
  #-----------------------------------------------------------------------------#

  try :
    silo = SiloOutput(mesh)
    silo.initialize("iaea3d.silo")
    silo.write_scalar_flux(state)
    silo.finalize()
  except :
    print "Silo error (not installed?)"

if __name__ == "__main__":
  Manager.initialize(sys.argv)
  run()
