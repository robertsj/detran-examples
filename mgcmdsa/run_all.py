"""
Archival CMDSA/MGDSA fixed-source comparisons for the C5G7 assembly.

The older branch for this example predated the Python user API and assumed
Python 2-era bindings.  This version keeps the same experiment structure while
compiling user-facing solver objects down to the legacy InputDB backend.
"""

from __future__ import print_function

import argparse
from dataclasses import dataclass
from pathlib import Path
import pickle
import sys
import time

import numpy as np
from detran import (
    ConstantSource,
    Fixed2D,
    IterationSolver,
    LinearSolver,
    Manager,
    Quadrature,
    ReactionRates,
    SolverDatabase,
    fixed_source_input,
)


SCRIPT_DIR = Path(__file__).resolve().parent
C5G7_DIR = SCRIPT_DIR.parent / "c5g7"
if str(C5G7_DIR) not in sys.path:
    sys.path.insert(0, str(C5G7_DIR))

from assemblies_c5g7 import get_assemblies
from material_c5g7 import get_materials


LEFT = 1
RIGHT = 2
GROUPS = 7
DEFAULT_QUADRATURE = 8
SMOKE_QUADRATURE = 1
_MANAGER_INITIALIZED = False


@dataclass(frozen=True)
class Case:
    key: str
    label: str
    code: tuple


ARCHIVAL_CASES = (
    Case("mg_gmres", "GMRES", (0, 1, 0, 0)),
    Case("mg_gmres_mgdsa_ilu0", "GMRES+MGDSA(ILU0)", (0, 1, 0, 1)),
    Case("mg_gmres_mgdsa_ilut", "GMRES+MGDSA(ILUT)", (0, 1, 0, 2)),
    Case("mg_gmres_cmdsa", "GMRES+CMDSA", (0, 1, 0, 3)),
    Case("mg_gmres_cmdsa_smooth", "GMRES+CMDSA+S", (0, 1, 0, 4)),
    Case("mg_gmres_tcmdsa_smooth", "GMRES+TC-CMDSA+S", (0, 1, 0, 5)),
    Case("cmfd", "CMFD", (1, 2, 0, 6)),
)

CONVERGENCE_CASES = (
    Case("mg_gmres", "GMRES", (0, 1, 0, 0)),
    Case("cmfd", "CMFD", (1, 2, 0, 6)),
)


def _vector(values):
    return np.asarray(list(values), dtype=float)


def _artifact(data_dir, name):
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir / name


def _backend_suffix(args):
    return "" if args.linear_backend == "builtin" else "_{}".format(args.linear_backend)


def _tol_label(tolerance):
    return "{:g}".format(float(tolerance))


def _quadrature(args, smoke=False):
    if args.quadrature is not None:
        return args.quadrature
    return SMOKE_QUADRATURE if smoke else DEFAULT_QUADRATURE


def _linear_database(
    *,
    rtol,
    atol=0.0,
    max_iters=5000,
    restart=30,
    monitor_level=0,
    pc_type=None,
    pc_side=LEFT,
    ilu_levels=None,
    ilu_thresh=None,
    backend="builtin",
):
    options = {}
    if backend == "builtin":
        if pc_type is not None:
            options["pc_type"] = pc_type
            options["pc_side"] = int(pc_side)
        if ilu_levels is not None:
            options["ilu_levels"] = int(ilu_levels)
        if ilu_thresh is not None:
            options["ilu_thresh"] = float(ilu_thresh)
        linear = LinearSolver.gmres(
            atol=atol,
            rtol=rtol,
            max_iters=max_iters,
            restart=restart,
            monitor_level=monitor_level,
            **options
        )
    elif backend == "petsc":
        options["linear_solver_gmres_restart"] = int(restart)
        if pc_type is not None:
            options["pc_type"] = "petsc_pc"
            options["pc_side"] = int(pc_side)
            if pc_type == "none":
                options["petsc_pc_type"] = "none"
            elif pc_type == "ilu0":
                options["petsc_pc_type"] = "ilu"
                options["petsc_pc_factor_levels"] = 0
            elif pc_type == "ilut":
                options["petsc_pc_type"] = "ilu"
                options["petsc_pc_factor_levels"] = int(ilu_levels or 2)
            else:
                raise ValueError("unsupported PETSc preconditioner: {}".format(pc_type))
        linear = LinearSolver.petsc(
            atol=atol,
            rtol=rtol,
            max_iters=max_iters,
            monitor_level=monitor_level,
            **options
        )
    else:
        raise ValueError("unsupported linear backend: {}".format(backend))
    return SolverDatabase(
        linear=linear
    )


def _preconditioner_database(pc_kind, rtol, backend="builtin"):
    if pc_kind == "none":
        return _linear_database(
            rtol=rtol,
            pc_type="none" if backend == "petsc" else None,
            backend=backend,
        )
    if pc_kind == "ilu0":
        return _linear_database(rtol=rtol, pc_type="ilu0", backend=backend)
    if pc_kind == "ilut":
        return _linear_database(
            rtol=rtol,
            pc_type="ilut",
            ilu_levels=2,
            ilu_thresh=0.0,
            backend=backend,
        )
    raise ValueError("unsupported preconditioner: {}".format(pc_kind))


def _inner_solver(code, tolerance):
    if code == 0:
        return IterationSolver.source_iteration(
            max_iters=20,
            tolerance=tolerance,
            print_level=0,
        )
    if code == 1:
        return IterationSolver.source_iteration(
            max_iters=1,
            tolerance=1.0e-16,
            print_level=0,
        )
    if code == 2:
        return IterationSolver.gmres(
            max_iters=20,
            tolerance=tolerance,
            print_level=0,
        )
    raise ValueError("invalid inner solver code: {}".format(code))


def _outer_solver(code, tolerance):
    if code == 0:
        return IterationSolver.gauss_seidel(
            max_iters=100000,
            tolerance=tolerance,
            print_level=2,
            print_interval=1,
        )
    if code == 1:
        return IterationSolver.gmres(
            max_iters=100000,
            tolerance=tolerance,
            print_level=2,
            print_interval=1,
        )
    if code == 2:
        return IterationSolver("CMFD", tolerance, 100000, 2, 1)
    raise ValueError("invalid outer solver code: {}".format(code))


def _apply_inner_pc(code, overrides):
    if code == 0:
        overrides["inner_pc_type"] = "none"
        return "none"
    if code == 1:
        overrides["inner_pc_type"] = "DSA"
        return "ilu0"
    if code == 2:
        overrides["inner_pc_type"] = "DSA"
        return "ilut"
    raise ValueError("invalid inner preconditioner code: {}".format(code))


def _apply_outer_pc(code, overrides):
    if code == 0:
        overrides["outer_pc_type"] = "none"
        return "none"
    if code == 1:
        overrides["outer_pc_type"] = "mgdsa"
        return "ilu0"
    if code == 2:
        overrides["outer_pc_type"] = "mgdsa"
        return "ilut"
    if code in (3, 4, 5):
        overrides["outer_pc_type"] = "mgcmdsa" if code in (3, 4) else "mgtcdsa"
        overrides["mgpc_coarse_mesh_level"] = 7
        overrides["mgpc_condensation_option"] = 1
        overrides["mgpc_cmdsa_use_smoothing"] = int(code in (4, 5))
        if code in (4, 5):
            overrides["mgpc_cmdsa_smoothing_iters"] = 3
            overrides["mgpc_cmdsa_smoothing_relax"] = 0.7
        if code == 5:
            overrides["mgpc_tcdsa_number_coarse_corrections"] = 3
            overrides["mgpc_tcdsa_use_fine_correction"] = 1
            overrides["mgpc_number_azimuth_octant"] = 2
            overrides["mgpc_number_polar_octant"] = 2
        return "ilut"
    if code == 6:
        overrides["mgpc_coarse_mesh_level"] = 7
        return "ilut"
    raise ValueError("invalid outer preconditioner code: {}".format(code))


def get_input(
    case,
    tolerance,
    quadrature_order=DEFAULT_QUADRATURE,
    pc_tolerance=0.01,
    linear_backend="builtin",
):
    inner_code, outer_code, inner_pc_code, outer_pc_code = tuple(case.code)
    overrides = {
        "store_angular_flux": 0,
        "store_current": int(outer_code == 2),
        "cmfd_relaxation": 1.0,
        "mgdsa_disable_fission": 0,
        "inner_pc_side": RIGHT,
        "outer_pc_side": RIGHT,
        "outer_krylov_group_cutoff": 0,
    }
    inner_pc = _apply_inner_pc(inner_pc_code, overrides)
    outer_pc = _apply_outer_pc(outer_pc_code, overrides)

    inner_solver_db = _linear_database(
        rtol=tolerance,
        max_iters=20 if inner_code == 2 else 50000,
        monitor_level=0,
        backend=linear_backend,
    )
    outer_solver_db = _linear_database(
        rtol=tolerance,
        max_iters=50000,
        monitor_level=2,
        backend=linear_backend,
    )
    inner_pc_db = _preconditioner_database(inner_pc, pc_tolerance, linear_backend)
    outer_pc_db = _preconditioner_database(
        outer_pc,
        1.0e-12 if outer_pc_code == 6 else pc_tolerance,
        linear_backend,
    )

    return fixed_source_input(
        groups=GROUPS,
        dimension=2,
        equation="dd",
        boundary="vacuum",
        inner=_inner_solver(inner_code, tolerance),
        outer=_outer_solver(outer_code, tolerance),
        quadrature=Quadrature(type="u-dgl", polar=quadrature_order, azimuth=quadrature_order),
        solver_databases={
            "inner_solver_db": inner_solver_db,
            "inner_pc_db": inner_pc_db,
            "outer_solver_db": outer_solver_db,
            "outer_pc_db": outer_pc_db,
        },
        input_overrides=overrides,
    )


def get_mesh():
    assemblies = get_assemblies(7, True)
    return assemblies[0].mesh()


def get_materials_with_diffusion():
    material = get_materials()
    material.compute_diff_coef()
    return material


def get_solver(
    case,
    tolerance,
    quadrature_order=DEFAULT_QUADRATURE,
    pc_tolerance=0.01,
    linear_backend="builtin",
):
    inp = get_input(case, tolerance, quadrature_order, pc_tolerance, linear_backend)
    material = get_materials_with_diffusion()
    mesh = get_mesh()
    solver = Fixed2D(inp, material, mesh, True)
    solver.setup()
    source = ConstantSource(GROUPS, solver.mesh(), 1.0, solver.quadrature())
    solver.set_source(source)
    solver.set_solver()
    return solver, material, mesh


def run_reference(args):
    quadrature_order = _quadrature(args)
    solver, material, mesh = get_solver(
        ARCHIVAL_CASES[0],
        1.0e-12,
        quadrature_order,
        args.pc_tolerance,
        args.linear_backend,
    )
    start = time.time()
    solver.solve()
    print("reference elapsed = {:.3f} s".format(time.time() - start))

    state = solver.state()
    rates = ReactionRates(material, mesh, state)
    reference_data = {
        "fission": _vector(rates.region_power("PINS", 1.0)),
        "quadrature": quadrature_order,
    }
    for group in range(GROUPS):
        reference_data["phi{}".format(group)] = _vector(state.phi(group))

    path = _artifact(args.data_dir, "reference_data{}.p".format(_backend_suffix(args)))
    with path.open("wb") as stream:
        pickle.dump(reference_data, stream, protocol=pickle.HIGHEST_PROTOCOL)
    print("wrote {}".format(path))
    return reference_data


def load_reference(args):
    path = _artifact(args.data_dir, "reference_data{}.p".format(_backend_suffix(args)))
    if not path.exists():
        raise FileNotFoundError(
            "missing {}; run `python3 run_all.py reference` first".format(path)
        )
    with path.open("rb") as stream:
        return pickle.load(stream)


def flux_norm(state, material, mesh, reference_data):
    norm = 0.0
    norm_reference = 0.0
    for group in range(GROUPS):
        phi = _vector(state.phi(group))
        reference_phi = reference_data["phi{}".format(group)]
        norm += np.linalg.norm(phi - reference_phi) ** 2
        norm_reference += np.linalg.norm(reference_phi) ** 2
    norm_phi = np.sqrt(norm) / np.sqrt(norm_reference)

    mask = reference_data["fission"] > 0.0
    rates = ReactionRates(material, mesh, state)
    fission = _vector(rates.region_power("PINS", 1.0))[mask]
    reference_fission = reference_data["fission"][mask]
    norm_fission = np.max(np.abs(fission - reference_fission) / reference_fission)
    return norm_phi, norm_fission


def flux_norm_vs_pin_power(args):
    quadrature_order = _quadrature(args)
    reference_data = load_reference(args)
    rows = []
    tolerance = 1.0
    for _ in range(10):
        tolerance *= 0.5
        solver, material, mesh = get_solver(
            ARCHIVAL_CASES[0],
            tolerance,
            quadrature_order,
            args.pc_tolerance,
            args.linear_backend,
        )
        solver.solve()
        norm_phi, norm_fission = flux_norm(solver.state(), material, mesh, reference_data)
        rows.append((tolerance, norm_phi, norm_fission))
        print("{:12.5e} {:12.5e} {:12.5e}".format(tolerance, norm_phi, norm_fission))

    path = _artifact(args.data_dir, "tol_data{}.txt".format(_backend_suffix(args)))
    np.savetxt(path, np.asarray(rows))
    print("wrote {}".format(path))


def run_convergence(args):
    quadrature_order = _quadrature(args)
    reference_data = load_reference(args)
    tolerances = {}
    sweeps = {}

    for case in CONVERGENCE_CASES:
        tolerance0 = 0.01
        tolerance1 = 0.001
        solver, material, mesh = get_solver(
            case,
            tolerance0,
            quadrature_order,
            args.pc_tolerance,
            args.linear_backend,
        )
        solver.solve()
        residual0, _ = flux_norm(solver.state(), material, mesh, reference_data)

        solver, material, mesh = get_solver(
            case,
            tolerance1,
            quadrature_order,
            args.pc_tolerance,
            args.linear_backend,
        )
        solver.solve()
        residual1, _ = flux_norm(solver.state(), material, mesh, reference_data)

        converged = False
        for iteration in range(50):
            if residual1 == residual0:
                tolerance = 1.0e-14
            else:
                slope = (tolerance1 - tolerance0) / (residual1 - residual0)
                tolerance = slope * (args.tolerance - residual1) + tolerance1
            tolerance = max(tolerance, 0.5 * tolerance1)

            tolerance0, residual0 = tolerance1, residual1
            tolerance1 = tolerance
            solver, material, mesh = get_solver(
                case,
                tolerance1,
                quadrature_order,
                args.pc_tolerance,
                args.linear_backend,
            )
            solver.solve()
            residual1, _ = flux_norm(solver.state(), material, mesh, reference_data)
            print(
                "{} iter {:02d}: input tolerance {:12.5e}, error {:12.5e}".format(
                    case.key,
                    iteration,
                    tolerance1,
                    residual1,
                )
            )
            if tolerance0 == tolerance1 or residual1 < args.tolerance:
                tolerances[case.code[0], case.code[1]] = tolerance1
                sweeps[case.key] = solver.number_sweeps()
                converged = True
                break
        if not converged:
            raise RuntimeError("could not calibrate {}".format(case.key))

    data = {"tolerances": tolerances, "sweeps": sweeps}
    path = _artifact(args.data_dir, "convergence_data{}.p".format(_backend_suffix(args)))
    with path.open("wb") as stream:
        pickle.dump(data, stream, protocol=pickle.HIGHEST_PROTOCOL)
    print("wrote {}".format(path))
    return data


def load_convergence(args):
    suffix = "" if args.calibration_source == "builtin" else _backend_suffix(args)
    path = _artifact(args.data_dir, "convergence_data{}.p".format(suffix))
    if not path.exists():
        raise FileNotFoundError(
            "missing {}; run `python3 run_all.py calibrate` first".format(path)
        )
    with path.open("rb") as stream:
        return pickle.load(stream)


def run_cases(args):
    quadrature_order = _quadrature(args)
    convergence = load_convergence(args)
    rows = []
    for case in ARCHIVAL_CASES:
        solver_tolerance = convergence["tolerances"][case.code[0], case.code[1]]
        print("{}: tolerance {}".format(case.key, solver_tolerance))
        solver, _, _ = get_solver(
            case,
            solver_tolerance,
            quadrature_order,
            args.pc_tolerance,
            args.linear_backend,
        )
        start = time.time()
        solver.solve()
        elapsed = time.time() - start
        sweeps = solver.number_sweeps()
        rows.append(
            {
                "case": case.key,
                "backend": args.linear_backend,
                "time": elapsed,
                "sweeps": sweeps,
            }
        )
        print("{:26s} {:10.3f} s {:8d} sweeps".format(case.label, elapsed, sweeps))

    path = _artifact(
        args.data_dir,
        "Xdata_{}_{}{}.p".format(
            _tol_label(args.tolerance),
            quadrature_order,
            _backend_suffix(args),
        ),
    )
    with path.open("wb") as stream:
        pickle.dump(rows, stream, protocol=pickle.HIGHEST_PROTOCOL)
    print("wrote {}".format(path))
    return rows


def load_case_data(args):
    quadrature_order = _quadrature(args)
    path = _artifact(
        args.data_dir,
        "Xdata_{}_{}{}.p".format(
            _tol_label(args.tolerance),
            quadrature_order,
            _backend_suffix(args),
        ),
    )
    if not path.exists():
        raise FileNotFoundError(
            "missing {}; run `python3 run_all.py run` first".format(path)
        )
    with path.open("rb") as stream:
        return pickle.load(stream)


def print_table(args):
    rows = load_case_data(args)
    by_key = {row["case"]: row for row in rows}

    out = []
    out.append("\\begin{table}[ht]")
    out.append("\\begin{center}")
    out.append("\\begin{tabular*}{0.7\\textwidth}{@{\\extracolsep{\\fill}} rccc}")
    out.append("\\toprule")
    out.append("                        & sweeps & sweeps/groups & time (s) \\\\")
    out.append("\\midrule")
    for case in ARCHIVAL_CASES:
        row = by_key[case.key]
        sweeps = int(row["sweeps"])
        out.append(
            "{} & {:4d} & {:4d} & {:8.1f} \\\\".format(
                case.label,
                sweeps,
                sweeps // GROUPS,
                float(row["time"]),
            )
        )
    out.append("\\bottomrule")
    out.append("\\end{tabular*}")
    out.append("\\end{center}")
    out.append(
        "\\caption{{Sweeps and timings for a C5G7 assembly with $\\tau = ${:8.4e} using {} linear solves}}".format(
            args.tolerance,
            args.linear_backend,
        )
    )
    out.append("\\end{table}")
    print("\n".join(out))


def smoke(args):
    quadrature_order = _quadrature(args, smoke=True)
    case = ARCHIVAL_CASES[3]
    solver, _, _ = get_solver(
        case,
        tolerance=0.5,
        quadrature_order=quadrature_order,
        pc_tolerance=args.pc_tolerance,
        linear_backend=args.linear_backend,
    )
    start = time.time()
    solver.solve()
    elapsed = time.time() - start
    print(
        "smoke {}: {:.3f} s, {} sweeps".format(
            case.label,
            elapsed,
            solver.number_sweeps(),
        )
    )


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        nargs="?",
        default="smoke",
        choices=("smoke", "reference", "tolerance-sweep", "calibrate", "run", "table", "all"),
    )
    parser.add_argument("--quadrature", type=int)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument("--pc-tolerance", type=float, default=0.01)
    parser.add_argument(
        "--linear-backend",
        choices=("builtin", "petsc"),
        default="builtin",
        help="linear solver backend for GMRES and diffusion preconditioner solves",
    )
    parser.add_argument(
        "--calibration-source",
        choices=("same", "builtin"),
        default="same",
        help="source for calibrated run tolerances",
    )
    parser.add_argument("--data-dir", default=str(SCRIPT_DIR))
    return parser.parse_args(argv)


def initialize_backend(args):
    global _MANAGER_INITIALIZED
    if args.linear_backend == "petsc" and not _MANAGER_INITIALIZED:
        Manager.initialize(sys.argv)
        _MANAGER_INITIALIZED = True


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    initialize_backend(args)
    if args.action == "smoke":
        smoke(args)
    elif args.action == "reference":
        run_reference(args)
    elif args.action == "tolerance-sweep":
        flux_norm_vs_pin_power(args)
    elif args.action == "calibrate":
        run_convergence(args)
    elif args.action == "run":
        run_cases(args)
    elif args.action == "table":
        print_table(args)
    elif args.action == "all":
        run_reference(args)
        if args.calibration_source == "same":
            run_convergence(args)
        run_cases(args)
        print_table(args)
    else:
        raise ValueError("unsupported action {}".format(args.action))


if __name__ == "__main__":
    main()
