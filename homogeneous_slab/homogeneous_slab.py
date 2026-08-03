"""
Here, a homogeneous, one-group slab is defined and analyzed using
a few different techniques.
"""

from detran import *
import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
import time

def get_solver_databases():
    solver_db = SolverDatabase(
        linear=LinearSolver.gmres(
            atol=1e-9,
            rtol=1e-9,
            max_iters=5000,
            restart=30,
            monitor_level=0,
            pc_type="ilu0",
            pc_side=1,
        )
    )
    return {
        "inner_solver_db": solver_db,
        "inner_pc_db": solver_db,
        "outer_solver_db": solver_db,
        "eigen_solver_db": solver_db,
    }

def get_input():
    """ fills default input parameters.
    """

    return fixed_source_input(
        groups=1,
        dimension=1,
        equation="dd",
        boundary="vacuum",
        inner=IterationSolver.source_iteration(
            max_iters=10000,
            tolerance=1e-12,
            print_level=0,
            use_pc=True,
        ),
        quadrature=Quadrature.gauss_legendre(1),
        solver_databases=get_solver_databases(),
        # For one group this is mostly the sweep count; keep the legacy output.
        input_overrides={"outer_print_out": 1},
    )

def get_mesh():
    # Mesh
    mesh = Mesh1D([200], [0.0, 100.0], [0])
    return mesh

def get_material(c=0.5):
    mat = Material(1, 1)
    mat.set_sigma_t(0, 0, 1.0)
    mat.set_sigma_s(0, 0, 0, c)
    mat.finalize()
    return mat

def run():
    mat = get_material()
    mesh = get_mesh()
    problem = FixedSourceProblem(
        mesh=mesh,
        materials=mat,
        boundary="vacuum",
        inner=IterationSolver.source_iteration(
            max_iters=10000,
            tolerance=1e-12,
            print_level=0,
            use_pc=True,
        ),
        quadrature=Quadrature.gauss_legendre(1),
        solver_databases=get_solver_databases(),
        input_overrides={"outer_print_out": 1},
    )
    solver = problem.create_raw_solver()
    solver.setup()
    quad = solver.quadrature()
    q_e = ConstantSource(1, mesh, 1.0, quad)
    solver.set_source(q_e)
    solver.set_solver()
    start = time.time()
    solver.solve()
    state = solver.state()
    elapsed = (time.time() - start)
    print(elapsed, " seconds")

    phi = np.asarray(state.phi(0))
    if plt is not None:
        edges = np.concatenate(
            ([0.0], np.cumsum([mesh.dx(i) for i in range(mesh.number_cells_x())]))
        )
        plt.stairs(phi, edges)
        plt.xlabel("x")
        plt.ylabel("scalar flux")
        plt.show()


if __name__ == "__main__":
  run()
