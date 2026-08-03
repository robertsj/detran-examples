"""
Consolidate kinetics input generation through the pydetran user API.
"""

from detran import (
    EigenSolver,
    IterationSolver,
    LinearSolver,
    Quadrature,
    SolverDatabase,
    TimeStepper,
    time_dependent_input,
)


def _case_settings(case):
    if case == "diffusion":
        return "diffusion", {"ts_discrete": 0, "store_angular_flux": 0}
    if case == "transport_full":
        return "dd", {"ts_discrete": 1, "store_angular_flux": 1}
    if case == "transport_iso":
        return "dd", {"ts_discrete": 0, "store_angular_flux": 0}
    raise ValueError("case must be diffusion, transport_full, or transport_iso")


def get_input(
    case="diffusion",
    *,
    adjoint=0,
    time_stepper=None,
    input_overrides=None,
):
    equation, case_overrides = _case_settings(case)
    solver_db = SolverDatabase(
        linear=LinearSolver.gmres(
            atol=0.0,
            rtol=1e-13,
            max_iters=1000,
            restart=50,
            monitor_level=0,
        ),
        eigen_algorithm="power",
        eigen_tolerance=1.0e-14,
        eigen_max_iters=1000,
        eigen_monitor_level=2,
    )
    overrides = {
        "adjoint": adjoint,
        "outer_krylov_group_cutoff": 0,
        "outer_pc_type": "none",
        "mgpc_coarse_mesh_level": 5,
        "mgpc_condensation_option": 1,
        "mgpc_cmdsa_use_smoothing": 1,
        "mgpc_cmdsa_smoothing_iters": 2,
        "mgpc_cmdsa_smoothing_relax": 0.8,
    }
    overrides.update(case_overrides)
    overrides.update(dict(input_overrides or {}))

    return time_dependent_input(
        groups=2,
        dimension=1,
        equation=equation,
        boundary="vacuum",
        solver=EigenSolver.power_iteration(
            tolerance=1e-14,
            max_iters=1000,
        ),
        inner=IterationSolver.gmres(
            tolerance=1e-12,
            max_iters=1000,
            print_level=0,
        ),
        outer=IterationSolver.gmres(
            tolerance=1e-12,
            max_iters=1000,
            print_level=0,
        ),
        quadrature=Quadrature(polar=8),
        time_stepper=time_stepper or TimeStepper.bdf1(
            step_size=0.05,
            final_time=60.0,
            max_steps=1000000,
            max_iters=1,
            tolerance=1.0e-5,
        ),
        solver_databases={
            "inner_solver_db": solver_db,
            "inner_pc_db": solver_db,
            "outer_solver_db": solver_db,
            "eigen_solver_db": solver_db,
        },
        input_overrides=overrides,
    )


def get_base_input():
    return get_input()


if __name__ == "__main__":
    inp = get_base_input()
    inp.display()
