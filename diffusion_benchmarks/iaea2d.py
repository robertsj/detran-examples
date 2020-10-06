"""
2-D IAEA benchmark

Reported keff ~ 1.02959

Reference:

Argonne Code Center. Benchmark  Problem  Book. Technical Report ANL-7416
Supplement 2, Argonne National Laboratory, 1977
"""

from detran import *
from diffusion_input import get_input
import time

def get_material():
    """
    Because this problem has the removal cross section specified, but
    not any within-group scattering, we'll simply put the removal values into
    the total, since removal is total minus within group.  Note also that the
    removal values are adjusted for axial buckling.
    """

    mat = Material.Create(4, 2, "IAEA-2D")
    # Material 0
    mat.set_sigma_t(0, vec_dbl([0.030120, 0.080032]))
    mat.set_sigma_s(0, 1, 0, 0.02)
    mat.set_diff_coef(0, vec_dbl([1.5, 0.4]))
    mat.set_sigma_f(0, 1, 0.135)
    mat.set_chi(0, 0, 1.0)
    # Material 1
    mat.set_sigma_t(1, vec_dbl([0.030120, 0.085032]))
    mat.set_sigma_s(1, 1, 0, 0.02)
    mat.set_diff_coef(1, vec_dbl([1.5, 0.4]))
    mat.set_sigma_f(1, 1, 0.135)
    mat.set_chi(1, 0, 1.0)
    # Material 2
    mat.set_sigma_t(2, vec_dbl([0.030120, 0.130032]))
    mat.set_sigma_s(2, 1, 0, 0.02)
    mat.set_diff_coef(2, vec_dbl([1.5, 0.4]))
    mat.set_sigma_f(2, 1, 0.135)
    mat.set_chi(2, 0, 1.0)
    ## Material 3
    mat.set_sigma_t(3, vec_dbl([0.040160, 0.010024]))
    mat.set_sigma_s(3, 1, 0, 0.04)
    mat.set_diff_coef(3, vec_dbl([2.0, 0.3]))
    mat.finalize()
    return mat

def get_mesh(num_div=1):
    """ Create mesh with num_div cells in x and y per quarter assembly.
    """
    fm = [num_div]+[2*num_div]*8
    cm = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0]
    mt = [2, 1, 1, 1, 2, 1, 1, 0, 3,
          1, 1, 1, 1, 1, 1, 1, 0, 3,
          1, 1, 1, 1, 1, 1, 0, 0, 3,
          1, 1, 1, 1, 1, 1, 0, 3, 3,
          2, 1, 1, 1, 2, 0, 0, 3, 3,
          1, 1, 1, 1, 0, 0, 3, 3, 3,
          1, 1, 0, 0, 0, 3, 3, 3, 3,
          0, 0, 0, 3, 3, 3, 3, 3, 3,
          3, 3, 3, 3, 3, 3, 3, 3, 3]
    mesh = Mesh2D.Create(fm, fm, cm, cm, mt)
    return mesh

def run():

    # Model
    inp = get_input()
    mat = get_material()
    mesh = get_mesh(5)

    # Solver
    solver = Eigen2D(inp, mat, mesh)
    t = time.time()
    solver.solve()
    state = solver.state()
    print "elapsed = ", time.time()-t

    return state.eigenvalue()

if __name__ == "__main__":
  Manager.initialize(sys.argv)
  run()
