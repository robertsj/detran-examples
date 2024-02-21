"""
2-D Biblis benchmark.

reported keff ~ 1.02513

Reference:


"""

from detran import *
from diffusion_input import get_input
import time

def get_material():
    """
    Note, because this problem has the removal cross section specified, but
    not any within-group scattering, we'll simply put the removal values into
    the total, since removal is total minus within group.
    """

    mat = Material(8, 2, "biblis")
    # Material 0
    mat.set_sigma_t(0, vec_dbl([0.0272582, 0.0750580]))
    mat.set_sigma_s(0, 1, 0, 0.017754)
    mat.set_diff_coef(0, vec_dbl([1.4360, 0.3635]))
    mat.set_sigma_f(0, vec_dbl([0.0058708, 0.0960670]))
    mat.set_chi(0, 0, 1.0)
    # Material 1
    mat.set_sigma_t(1, vec_dbl([0.0272995, 0.0784360]))
    mat.set_sigma_s(1, 1, 0, 0.017621)
    mat.set_diff_coef(1, vec_dbl([1.4366, 0.3636]))
    mat.set_sigma_f(1, vec_dbl([0.0061908, 0.1035800]))
    mat.set_chi(1, 0, 1.0)
    # Material 2
    mat.set_sigma_t(2, vec_dbl([0.0257622, 0.0715960]))
    mat.set_sigma_s(2, 1, 0, 0.023106)
    mat.set_diff_coef(2, vec_dbl([1.3200, 0.2772]))
    # Material 3
    mat.set_sigma_t(3, vec_dbl([0.0274640, 0.0914080]))
    mat.set_sigma_s(3, 1, 0, 0.017101)
    mat.set_diff_coef(3, vec_dbl([1.4389, 0.3638]))
    mat.set_sigma_f(3, vec_dbl([0.0074527, 0.1323600]))
    mat.set_chi(3, 0, 1.0)
    # Material 4
    mat.set_sigma_t(4, vec_dbl([0.0272930, 0.0848280]))
    mat.set_sigma_s(4, 1, 0, 0.017290)
    mat.set_diff_coef(4, vec_dbl([1.4381, 0.3665]))
    mat.set_sigma_f(4, vec_dbl([0.0061908, 0.1035800]))
    mat.set_chi(4, 0, 1.0)
    # Material 5
    mat.set_sigma_t(5, vec_dbl([0.0273240, 0.0873140]))
    mat.set_sigma_s(5, 1, 0, 0.017192)
    mat.set_diff_coef(5, vec_dbl([1.4385, 0.3665]))
    mat.set_sigma_f(5, vec_dbl([0.0064285, 0.1091100]))
    mat.set_chi(5, 0, 1.0)
    # Material 6
    mat.set_sigma_t(6, vec_dbl([0.0272900, 0.0880240]))
    mat.set_sigma_s(6, 1, 0, 0.017125)
    mat.set_diff_coef(6, vec_dbl([1.4389, 0.3679]))
    mat.set_sigma_f(6, vec_dbl([0.0061908, 0.1035800]))
    mat.set_chi(6, 0, 1.0)
    # Material 7
    mat.set_sigma_t(7, vec_dbl([0.0273210, 0.0905100]))
    mat.set_sigma_s(7, 1, 0, 0.017027)
    mat.set_diff_coef(7, vec_dbl([1.4393, 0.3680]))
    mat.set_sigma_f(7, vec_dbl([0.0064285, 0.1091100]))
    mat.set_chi(7, 0, 1.0)
    mat.finalize()
    return mat

def get_mesh(num_div=1):
    """ Create mesh with num_div cells in x and y per quarter assembly.
    """
    cm = [0.0, 11.5613]
    for i in range(2, 10) :
        cm.append(cm[-1] + 23.1226)
    fm = [num_div] + [2*num_div]*8
    mt = [0, 7, 1, 5, 0, 6, 0, 3, 2,
          7, 0, 7, 1, 7, 0, 0, 3, 2,
          1, 7, 0, 7, 1, 6, 0, 3, 2,
          5, 1, 7, 1, 7, 0, 7, 3, 2,
          0, 7, 1, 7, 1, 4, 3, 2, 2,
          6, 0, 6, 0, 4, 3, 3, 2, 2,
          0, 0, 0, 7, 3, 3, 2, 2, 2,
          3, 3, 3, 3, 2, 2, 2, 2, 2,
          2, 2, 2, 2, 2, 2, 2, 2, 2]
    mesh = Mesh2D(fm, fm, cm, cm, mt)
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
    print("elapsed = ", time.time()-t)
    return state.eigenvalue()

if __name__ == "__main__":
  #Manager.initialize(sys.argv)
  run()
