"""
Base input database for diffusion problems.
"""

from detran import (
    Boundary,
    DiffusionEigenvalue,
    EigenSolver,
    LinearSolver,
    diffusion_eigenvalue_input,
)

DEFAULT_BOUNDARY = Boundary(
    west="reflect",
    east="vacuum",
    south="reflect",
    north="vacuum",
    bottom="vacuum",
    top="vacuum",
)


def get_solver(solver=None):
    """Default solver parameters for diffusion benchmark eigenproblems."""

    if isinstance(solver, EigenSolver):
        return solver

    linear = LinearSolver.gmres(
        atol=1e-8,
        rtol=1e-8,
        max_iters=5000,
        restart=30,
        monitor_level=0,
    )

    if solver == "slepc":
        return EigenSolver.slepc(
            inner=linear,
            max_iters=1000,
            monitor_level=0,
        )

    return EigenSolver.arnoldi(
        algorithm="gd",
        inner=linear,
        max_iters=1000,
        monitor_level=0,
    )


def get_input(solver=None, dimension=2):
    """Compatibility helper returning the generated legacy InputDB."""

    return diffusion_eigenvalue_input(
        groups=2,
        dimension=dimension,
        boundary=DEFAULT_BOUNDARY,
        solver=get_solver(solver),
    )


def get_problem(mesh, material, solver=None, dimension=None):
    """Build a user-facing diffusion eigenvalue problem."""

    return DiffusionEigenvalue(
        mesh=mesh,
        materials=material,
        boundary=DEFAULT_BOUNDARY,
        solver=get_solver(solver),
        dimension=dimension,
    )




if __name__ == "__main__":
    inp = get_input()
    inp.display()
