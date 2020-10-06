# pyexamples/slab_reactor/slab_reactor.py
#
# This implements the test cases published in Scott Mosher's
# Ph.D. thesis, "A Variational Coarse Mesh Transport Method".
# All data and reference values are from Appendix A of that
# work.
#
#   Core     keff
#   =====   =========
#     1     1.258247
#     2     1.007066
#     3     0.805372
#
# These results are based on a 32-point Gauss-Legendre,
# 2 meshes per water region, 4 meshes per fuel region.
#
# Note, the assembly kinf's were found using the same
# parameters to be
#
#   Assembly     kinf
#   ========   =========
#      A       1.329576914
#      B       1.298737436
#      C       0.681362819
#      D       0.191909997

from detran import *
import time

def get_materials() :
    # Two-group data from 1-d coarse mesh benchmarks (Mosher, Ilas etc.)

    # Create the Materials object.
    mat = Material.Create(4, 2, "slabreactor");

    # ---------------------------
    # Material 0: Water
    # ---------------------------

    # Total
    mat.set_sigma_t(0, 0, 0.1890);       # (obj, matid, g, value);
    mat.set_sigma_t(0, 1, 1.4633);

    # Fission
    mat.set_sigma_f(0, 0, 0.0);         # Note, default is zero
    mat.set_sigma_f(0, 1, 0.0);
    mat.set_chi(0, 0, 0.0);
    mat.set_chi(0, 1, 0.0);

    # Scattering
    mat.set_sigma_s(0, 0, 0, 0.1507);    # 1 <- 1
    mat.set_sigma_s(0, 0, 1, 0.0000);    # 1 <- 2
    mat.set_sigma_s(0, 1, 0, 0.0380);    # 2 <- 1
    mat.set_sigma_s(0, 1, 1, 1.4536);    # 2 <- 2

    # ---------------------------
    # Material 1: Fuel I
    # ---------------------------

    # Total
    mat.set_sigma_t(1, 0, 0.2263);       # (obj, matid, g, value);
    mat.set_sigma_t(1, 1, 1.0119);

    # Fission
    mat.set_sigma_f(1, 0, 0.0067);
    mat.set_sigma_f(1, 1, 0.1241);
    mat.set_chi(1, 0, 1.0);
    mat.set_chi(1, 1, 0.0);

    # Scattering
    mat.set_sigma_s(1, 0, 0, 0.2006);    # 1 <- 1
    mat.set_sigma_s(1, 0, 1, 0.0000);    # 1 <- 2
    mat.set_sigma_s(1, 1, 0, 0.0161);    # 2 <- 1
    mat.set_sigma_s(1, 1, 1, 0.9355);    # 2 <- 2

    # ---------------------------
    # Material 3: Fuel II
    # ---------------------------

    # Total
    mat.set_sigma_t(2, 0, 0.2252);       # (obj, matid, g, value);
    mat.set_sigma_t(2, 1, 0.9915);

    # Fission
    mat.set_sigma_f(2, 0, 0.0078);
    mat.set_sigma_f(2, 1, 0.1542);
    mat.set_chi(2, 0, 1.0);
    mat.set_chi(2, 1, 0.0);

    # Scattering
    mat.set_sigma_s(2, 0, 0, 0.1995);    # 1 <- 1
    mat.set_sigma_s(2, 0, 1, 0.0000);    # 1 <- 2
    mat.set_sigma_s(2, 1, 0, 0.0156);    # 2 <- 1
    mat.set_sigma_s(2, 1, 1, 0.9014);    # 2 <- 2

    # ---------------------------
    # Material 4: Fuel II + Gd
    # ---------------------------

    # Total
    mat.set_sigma_t(3, 0, 0.2173);       # (obj, matid, g, value);
    mat.set_sigma_t(3, 1, 1.0606);

    # Fission
    mat.set_sigma_f(3, 0, 0.0056);
    mat.set_sigma_f(3, 1, 0.0187);
    mat.set_chi(3, 0, 1.0);
    mat.set_chi(3, 1, 0.0);

    # Scattering
    mat.set_sigma_s(3, 0, 0, 0.1902);	   # 1 <- 1
    mat.set_sigma_s(3, 0, 1, 0.0000);    # 1 <- 2
    mat.set_sigma_s(3, 1, 0, 0.0136);    # 2 <- 1
    mat.set_sigma_s(3, 1, 1, 0.5733);    # 2 <- 2

    # ---------------------------
    # FINALIZE
    # ---------------------------

    mat.finalize();
    return mat

def get_mesh(geometry_id) :

  # Define the assembly coarse mesh edges.
  cm_assembly = [1.1580, 4.4790, 7.8000, 11.1210, 14.4420, 15.6000]

  # Define the base fine mesh count per coarse mesh.
  fm_assembly = [2, 4, 4, 4, 4, 2]

  # A core is composed of 7 adjacent assemblies.
  cm_core = [0.0]
  fm_core = []
  for i in range(0, 7) :
    for j in range(0, len(cm_assembly)) :
      cm_core.append(cm_assembly[j] + 15.6*float(i))
      fm_core.append(fm_assembly[j])

  # Add the 0.0 edge to assembly coarse mesh
  cm_assembly.insert(0, 0.0)

  # Define the coarse mesh material maps for each assembly type.
  # There are four of these assemblies.  These correspond, in order,
  # to types A, B, C, and D in the thesis.
  assem = [[ 0, 1, 2, 2, 1, 0 ], \
           [ 0, 1, 1, 1, 1, 0 ], \
           [ 0, 1, 3, 3, 1, 0 ], \
           [ 0, 3, 3, 3, 3, 0 ]]

  # Cores 0, 1 and 2
  core_0 = assem[0]+assem[1]+assem[0]+assem[1]+assem[0]+assem[1]+assem[0]
  core_1 = assem[0]+assem[2]+assem[0]+assem[2]+assem[0]+assem[2]+assem[0]
  core_2 = assem[0]+assem[3]+assem[0]+assem[3]+assem[0]+assem[3]+assem[0]

  # Create the 1D mesh
  if   geometry_id == "assembly0" :
    mesh = Mesh1D.Create(fm_assembly, cm_assembly, assem[0])
  elif geometry_id == "assembly1" :
    mesh = Mesh1D.Create(fm_assembly, cm_assembly, assem[1])
  elif geometry_id == "assembly2" :
    mesh = Mesh1D.Create(fm_assembly, cm_assembly, assem[2])
  elif geometry_id == "assembly3" :
    mesh = Mesh1D.Create(fm_assembly, cm_assembly, assem[3])
  elif geometry_id == "core0" :
    mesh = Mesh1D.Create(fm_core, cm_core, core_0)
  elif geometry_id == "core1" :
    mesh = Mesh1D.Create(fm_core, cm_core, core_1)
  elif geometry_id == "core2" :
    mesh = Mesh1D.Create(fm_core, cm_core, core_2)
  else :
    print "invalid geometry selected."
    exit()

  return mesh


def get_input():
    inp = InputDB.Create()
    inp.put_str("problem_type",               "eigenvalue")
    inp.put_int("number_groups",              2)
    inp.put_str("equation",                   "dd")
    inp.put_str("inner_solver",               "SI")
    inp.put_int("inner_max_iters",            1000)
    inp.put_dbl("inner_tolerance",            1e-7)
    inp.put_int("inner_print_level",          0)
    inp.put_str("outer_solver",               "GS")
    inp.put_int("outer_max_iters",            1000)
    inp.put_dbl("outer_tolerance",            1e-7)
    inp.put_int("outer_print_level",          0)
    inp.put_str("eigen_solver",               "PI")
    inp.put_int("eigen_max_iters",            200)
    inp.put_dbl("eigen_tolerance",            1e-7)
    inp.put_str("bc_west",                    "vacuum")
    inp.put_str("bc_east",                    "vacuum")
    inp.put_str("quad_type",                  "gl")
    inp.put_int("quad_number_polar_octant",   16)
    return inp

def run():
    inp = get_input()
    mat = get_materials()

    cases = ["core0", "core1", "core2",
             "assembly0", "assembly1", "assembly2", "assembly3"]
    keffs = []
    for case in cases:
        if "assem" in case:
            inp.put_str("bc_west", "reflect")
            inp.put_str("bc_east", "reflect")
        mesh = get_mesh(case)
        solver = Eigen1D(inp, mat, mesh)
        solver.solve()
        keffs.append(solver.state().eigenvalue())
    return cases, keffs

if __name__ == "__main__":
    Manager.initialize(sys.argv)
    run()
