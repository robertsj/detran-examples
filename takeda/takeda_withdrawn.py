# Takeda benchmark, model 1, inserted control (case 2)
# Reference keff ~ 0.9780
# 3-D NEUTRON TRANSPORT BENCHMARKS, NEACRP-L-330

from detran import *
import time

def get_input():
  solver_db = SolverDatabase(
    linear=LinearSolver.gmres(
      atol=1e-8,
      rtol=1e-8,
      max_iters=5000,
      restart=30,
      monitor_level=0,
    ),
    eigen_algorithm="gd",
    eigen_monitor_level=0,
  )
  return eigenvalue_input(
    groups=2,
    dimension=3,
    equation="dd",
    boundary=("reflect", "vacuum", "reflect", "vacuum", "reflect", "vacuum"),
    inner=IterationSolver.gmres(
      max_iters=1,
      tolerance=1e-6,
      print_level=0,
      print_interval=10,
      use_pc=False,
      pc_tolerance=1e-2,
      pc_ilu_levels=0,
      pc_type="none",
    ),
    outer=IterationSolver.gmres(
      max_iters=0,
      tolerance=1e-6,
      print_level=0,
      print_interval=2,
      upscatter_cutoff=0,
      use_pc=False,
    ),
    solver=EigenSolver.arnoldi(
      max_iters=1000,
      tolerance=1e-6,
      monitor_level=2,
      print_interval=1,
    ),
    quadrature=Quadrature(polar=2, azimuth=2),
    solver_databases={
      "inner_solver_db": solver_db,
      "outer_solver_db": solver_db,
      "eigen_solver_db": solver_db,
    },
  )

def run() :

  inp = get_input()

  # Mesh
  cm_x = [0.0, 15.0, 20.0, 25.0]
  cm_y = [0.0,  5.0, 15.0, 25.0]
  cm_z = [0.0, 15.0, 25.0]
  fm_x = [30, 10, 10]
  fm_y = [10, 20, 20]
  fm_z = [30, 20]
  mt   = [# Fueled region
          0, 2, 1,
          0, 1, 1,
          1, 1, 1,
          # Above fuel
          1, 2, 1,
          1, 1, 1,
          1, 1, 1]

  mesh = Mesh3D(fm_x, fm_y, fm_z, cm_x, cm_y, cm_z, mt)

  # Material
  mat = Material(3, 2)
  # 0. Core
  mat.set_sigma_t(0, 0,    2.23775e-01)
  mat.set_sigma_t(0, 1,    1.03864e-00)
  mat.set_sigma_s(0, 0, 0, 1.92423e-01)
  mat.set_sigma_s(0, 1, 0, 2.28253e-02)
  mat.set_sigma_s(0, 1, 1, 8.80439e-01)
  mat.set_sigma_f(0, 0,    9.09319e-03)
  mat.set_sigma_f(0, 1,    2.90183e-01)
  mat.set_chi(0, 0,        1.0)
  # 1. Reflector
  mat.set_sigma_t(1, 0,    2.50367e-01)
  mat.set_sigma_t(1, 1,    1.64482e-00)
  mat.set_sigma_s(1, 0, 0, 1.93446e-01)
  mat.set_sigma_s(1, 1, 0, 5.65042e-02)
  mat.set_sigma_s(1, 1, 1, 1.62452e-00)
  # 2. Control Rod (Widthdrawn = Void)
  mat.set_sigma_t(2, 0,    1.28407e-02)
  mat.set_sigma_t(2, 1,    1.20676e-02)
  mat.set_sigma_s(2, 0, 0, 1.27700e-02)
  mat.set_sigma_s(2, 1, 0, 2.40997e-05)
  mat.set_sigma_s(2, 1, 1, 1.07387e-02)
  #
  mat.finalize()

  try:
    io = IO_HDF5("takeda_withdrawn.h5")
    io.open()
    io.write(inp)
    io.write(mat)
    io.write(mesh)
    io.close()
  except:
    print("Error using HDF5---perhaps detran is not built with it.")

  # Solve
  start = time.time()
  solver = Eigen3D(inp, mat, mesh)
  solver.solve()
  print("elapsed = ", time.time() - start)
  state = solver.state()

  # Save fluxes
  try :
    silo = SiloOutput(mesh)
    silo.initialize("takeda_withdrawn.silo")
    silo.write_mesh_map("MATERIAL")
    silo.write_scalar_flux(state)
    silo.finalize()
  except :
    print("Error using Silo---perhaps detran is not built with it.")

  return state.eigenvalue()

if __name__ == "__main__":
  run()
