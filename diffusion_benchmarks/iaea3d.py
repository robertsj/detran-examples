"""
Solves the 3D IAEA PWR benchmark problem.  Note, because of modeling
limitations, the void regions beyond the normal reflector are simply
filled with more reflector.

reported keff ~ 1.029096

Reference:
  Benchmark Problem Book, ANL-7416, Suppl. 2, Argonne National
  Laboratory (1977)
"""
import sys
from detran import *
from diffusion_input import get_input
import time

def get_material():
    """
    All absorption cross sections are simply put into the total.
    """

    mat = Material(6, 2, "IAEA-3D")
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
    return mat

def get_mesh(num_div_xy=1, num_div_z=1):
    """ xy divs per 10 cm and z divs per 10 cm.
    """
    # This sets up for a 2cm mesh in all directions.
    # XY plane discretization
    cmH = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0]
    fmH = [num_div_xy]+[2*num_div_xy]*8
    # Axial discretization
    cmV = [0.0, 20.0, 280.0, 360.0, 380.0]
    fmV = [num_div_z*i for i in [2,26,8,2]]
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
    mesh = Mesh3D(fmH, fmH, fmV, cmH, cmH, cmV, mt)
    return mesh

def run():

    inp = get_input()
    inp.put_int("dimension", 3)
    mat = get_material()
    mesh = get_mesh(2, 2)

    # Solver
    solver = Eigen3D(inp, mat, mesh)
    t = time.time()
    solver.solve()
    state = solver.state()
    print("elapsed = ", time.time()-t)
    return state.eigenvalue()

if __name__ == "__main__":
  #Manager.initialize(sys.argv)
  run()
