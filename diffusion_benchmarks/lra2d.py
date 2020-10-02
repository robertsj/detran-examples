"""
2-D LRA benchmark.

reported keff ~ 0.99636

Reference:

"""

from detran import *
from diffusion_input import get_input
import time


def get_material():
    """
    Note, because this problem has the removal cross section specified, but
    not any within-group scattering, we'll simply put the removal values into
    the total, since removal is total minus within group.  Note also that the
    removal values are adjusted for axial buckling.
    """

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
    return mat

def get_mesh(num_div=1):
    """ Create mesh with num_div cells per side per assembly (15x15 cm)
    """
    cm = [0.0, 15.0, 75.0, 105.0, 120.0, 135.0, 165.0]
    fm = [num_div*i for i in [1, 4, 2, 1, 1, 2]]
    mt = [1, 0, 1, 2, 2, 4,
          0, 0, 0, 2, 2, 4,
          1, 0, 1, 2, 2, 4,
          2, 2, 2, 3, 4, 4,
          2, 2, 2, 4, 4, 4,
          4, 4, 4, 4, 4, 4]
    mesh = Mesh2D.Create(fm, fm, cm, cm, mt)
    return mesh

def run():

    # Model
    inp = get_input()
    mat = get_material()
    mesh = get_mesh()

    # Solver
    solver = Eigen2D(inp, mat, mesh)
    t = time.time()
    solver.solve()
    state = solver.state()
    print "elapsed = ", time.time()-t

if __name__ == "__main__":
  Manager.initialize(sys.argv)
  run()
