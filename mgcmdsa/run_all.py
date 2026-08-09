"""
Archival CMDSA/MGDSA fixed-source comparisons for the C5G7 assembly.

The older branch for this example predated the Python user API and assumed
Python 2-era bindings.  This version keeps the same experiment structure while
compiling user-facing solver objects down to the legacy InputDB backend.
"""

from __future__ import print_function

import argparse
import copy
from dataclasses import dataclass
import json
from pathlib import Path
import pickle
import sys
import time

import numpy as np
import detran as _detran
from detran import (
    ConstantSource,
    FixedBoundarySource,
    FixedSourceProblem,
    IterationSolver,
    LinearSolver,
    Manager,
    Material,
    Mesh2D,
    Quadrature,
    ReactionRates,
    SolverDatabase,
    fixed_source_input,
)


SCRIPT_DIR = Path(__file__).resolve().parent
PAPER_RESULTS = SCRIPT_DIR / "paper_results.json"
UPDATED_RESULTS_DIR = SCRIPT_DIR / "updated_results"
UPDATED_RESULTS_C5G7 = "c5g7_local_runs.json"
UPDATED_RESULTS_BWR_DIAGNOSTICS = "bwr_spectral_diagnostics_with_fission.json"
UPDATED_RESULTS_BWR_EIGS = "bwr_spectral_eigenvalues_with_fission.npz"
UPDATED_RESULTS_REPORT = "updated_results.pdf"
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
    paper_key: str = None


PAPER_C5G7_CASES = (
    Case("gs_si20", "GS, SI(20)", (0, 0, 0, 0), "gs_si20"),
    Case("gs_si1", "GS, SI(1)", (1, 0, 0, 0), "gs_si1"),
    Case("gs_gmres", "GS, GMRES", (2, 0, 0, 0), "gs_gmres"),
    Case("gs_dsa_lu", "GS, GMRES+DSA(LU)", (2, 0, 1, 0), "gs_dsa_lu"),
    Case("gs_dsa_ilu", "GS, GMRES+DSA(ILU)", (2, 0, 2, 0), "gs_dsa_ilu"),
    Case("gmres", "GMRES", (0, 1, 0, 0), "gmres"),
    Case("dsa_lu", "GMRES+DSA(LU)", (0, 1, 0, 1), "dsa_lu"),
    Case("dsa_ilu", "GMRES+DSA(ILU)", (0, 1, 0, 2), "dsa_ilu"),
    Case("cmdsa", "GMRES+CMDSA", (0, 1, 0, 3), "cmdsa"),
    Case("cmdsa_s", "GMRES+CMDSA+S", (0, 1, 0, 4), "cmdsa_s"),
    Case("tc1_cmdsa_s", "GMRES+TC-CMDSA+S", (0, 1, 0, 5), "tc1_cmdsa_s"),
)

ARCHIVAL_CASES = PAPER_C5G7_CASES

CONVERGENCE_CASES = (
    Case("gs_si20", "GS, SI(20)", (0, 0, 0, 0), "gs_si20"),
    Case("gs_si1", "GS, SI(1)", (1, 0, 0, 0), "gs_si1"),
    Case("gs_gmres", "GS, GMRES", (2, 0, 0, 0), "gs_gmres"),
    Case("gmres", "GMRES", (0, 1, 0, 0), "gmres"),
)

CURRENT_TO_PAPER_CASES = {case.key: case.paper_key for case in ARCHIVAL_CASES}

BWR_OPERATOR_KINDS = ("transport", "mgdsa", "cmdsa", "cmdsa-s", "tcmdsa")
BWR_REPORT_ORDER = ("A", "DSA", "CMDSA", "CMDSA-S", "TC-CMDSA")
BWR_EIGENVALUE_SOURCES = {
    "A": ("A", "eigvals_A_with_fission.npz"),
    "DSA": ("DSA", "eigvals_DSA_with_fission.npz"),
    "CMDSA": ("CMDSA", "eigvals_CMDSA_with_fission.npz"),
    "CMDSA-S": ("CMDSA_S", "eigvals_CMDSA-S_with_fission.npz"),
}


def case_by_key(key):
    for case in ARCHIVAL_CASES:
        if case.key == key:
            return case
    raise KeyError(key)


def _selected_cases(args, cases):
    if not args.cases:
        return tuple(cases)
    by_key = {case.key: case for case in cases}
    selected = []
    for key in args.cases:
        if key not in by_key:
            raise ValueError(
                "case {!r} is not valid for this action; choose from {}".format(
                    key,
                    ", ".join(case.key for case in cases),
                )
            )
        selected.append(by_key[key])
    return tuple(selected)


def _ordered_rows(rows):
    by_key = {row["case"]: row for row in rows}
    ordered = [by_key[case.key] for case in ARCHIVAL_CASES if case.key in by_key]
    extras = [row for row in rows if row["case"] not in CURRENT_TO_PAPER_CASES]
    return ordered + extras


def _vector(values):
    return np.asarray(list(values), dtype=float)


def _artifact(data_dir, name):
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir / name


def _reference_path(args, quadrature_order=None):
    if quadrature_order is None:
        quadrature_order = _quadrature(args)
    return _artifact(
        args.data_dir,
        "reference_data_q{}{}.p".format(
            int(quadrature_order),
            _backend_suffix(args),
        ),
    )


def _legacy_reference_path(args):
    return _artifact(args.data_dir, "reference_data{}.p".format(_backend_suffix(args)))


def _convergence_path(args, quadrature_order=None):
    if quadrature_order is None:
        quadrature_order = _quadrature(args)
    return _artifact(
        args.data_dir,
        "convergence_data_{}_q{}{}.p".format(
            _tol_label(args.tolerance),
            int(quadrature_order),
            _backend_suffix(args),
        ),
    )


def _legacy_convergence_path(args):
    suffix = "" if args.calibration_source == "builtin" else _backend_suffix(args)
    return _artifact(args.data_dir, "convergence_data{}.p".format(suffix))


def _case_data_path(args, quadrature_order=None):
    if quadrature_order is None:
        quadrature_order = _quadrature(args)
    return _artifact(
        args.data_dir,
        "Xdata_{}_{}{}.p".format(
            _tol_label(args.tolerance),
            quadrature_order,
            _backend_suffix(args),
        ),
    )


def _tol_sweep_path(args, quadrature_order=None):
    if quadrature_order is None:
        quadrature_order = _quadrature(args)
    return _artifact(
        args.data_dir,
        "tol_data_q{}{}.txt".format(
            int(quadrature_order),
            _backend_suffix(args),
        ),
    )


def _legacy_tol_sweep_path(args):
    return _artifact(args.data_dir, "tol_data{}.txt".format(_backend_suffix(args)))


def _backend_suffix(args):
    return "" if args.linear_backend == "builtin" else "_{}".format(args.linear_backend)


def _tol_label(tolerance):
    return "{:g}".format(float(tolerance))


def _paper_tol_key(tolerance):
    label = "{:.0e}".format(float(tolerance))
    label = label.replace("e-0", "e-").replace("e+0", "e")
    return label.replace("e+", "e")


def _quadrature(args, smoke=False):
    if args.quadrature is not None:
        return args.quadrature
    return SMOKE_QUADRATURE if smoke else DEFAULT_QUADRATURE


def load_paper_results(path=PAPER_RESULTS):
    with Path(path).open("r", encoding="utf-8") as stream:
        return json.load(stream)


def _paper_residual_seed(args, case):
    paper_key = CURRENT_TO_PAPER_CASES.get(case.key)
    if paper_key is None:
        return None
    tol_key = _paper_tol_key(args.tolerance)
    tables = load_paper_results(args.paper_results)["c5g7_fixed_source"]["tables"]
    for row in tables.get(tol_key, []):
        if row["key"] == paper_key:
            return row["residual"]
    return None


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
            if pc_type == "lu":
                pc_type = "ilut"
                ilu_levels = 10000 if ilu_levels is None else ilu_levels
                ilu_thresh = 0.0 if ilu_thresh is None else ilu_thresh
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
            elif pc_type == "lu":
                options["petsc_pc_type"] = "lu"
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
    if pc_kind == "lu":
        return _linear_database(
            rtol=rtol,
            pc_type="lu",
            ilu_levels=10000,
            ilu_thresh=0.0,
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
        return "lu"
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
        return "lu"
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


def case_components(
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
    solver_databases = {
        "inner_solver_db": inner_solver_db,
        "inner_pc_db": inner_pc_db,
        "outer_solver_db": outer_solver_db,
        "outer_pc_db": outer_pc_db,
    }
    quadrature = Quadrature(
        type="u-dgl",
        polar=quadrature_order,
        azimuth=quadrature_order,
    )
    return {
        "inner": _inner_solver(inner_code, tolerance),
        "outer": _outer_solver(outer_code, tolerance),
        "quadrature": quadrature,
        "solver_databases": solver_databases,
        "input_overrides": overrides,
    }


def get_input(
    case,
    tolerance,
    quadrature_order=DEFAULT_QUADRATURE,
    pc_tolerance=0.01,
    linear_backend="builtin",
):
    components = case_components(
        case,
        tolerance,
        quadrature_order,
        pc_tolerance,
        linear_backend,
    )

    return fixed_source_input(
        groups=GROUPS,
        dimension=2,
        equation="dd",
        boundary={"west": "fixed", "east": "vacuum", "south": "vacuum", "north": "vacuum"},
        inner=components["inner"],
        outer=components["outer"],
        quadrature=components["quadrature"],
        solver_databases=components["solver_databases"],
        input_overrides=components["input_overrides"],
    )


def get_mesh():
    assemblies = get_assemblies(7, True)
    return assemblies[0].mesh()


def get_materials_with_diffusion():
    material = get_materials()
    material.compute_diff_coef()
    return material


def c5g7_boundary_source():
    return FixedBoundarySource(side="west", group=0, value=1.0)


def get_problem(
    case,
    tolerance,
    quadrature_order=DEFAULT_QUADRATURE,
    pc_tolerance=0.01,
    linear_backend="builtin",
):
    components = case_components(
        case,
        tolerance,
        quadrature_order,
        pc_tolerance,
        linear_backend,
    )
    material = get_materials_with_diffusion()
    mesh = get_mesh()
    problem = FixedSourceProblem(
        mesh=mesh,
        materials=material,
        boundary={"west": "fixed", "east": "vacuum", "south": "vacuum", "north": "vacuum"},
        equation="dd",
        inner=components["inner"],
        outer=components["outer"],
        quadrature=components["quadrature"],
        solver_databases=components["solver_databases"],
        input_overrides=components["input_overrides"],
        multiply=True,
        fixed_boundary_sources=[c5g7_boundary_source()],
    )
    return problem, material, mesh


def get_solver(
    case,
    tolerance,
    quadrature_order=DEFAULT_QUADRATURE,
    pc_tolerance=0.01,
    linear_backend="builtin",
):
    problem, material, mesh = get_problem(
        case,
        tolerance,
        quadrature_order,
        pc_tolerance,
        linear_backend,
    )
    return problem.prepare_raw_solver(), material, mesh


def run_reference(args):
    quadrature_order = _quadrature(args)
    solver, material, mesh = get_solver(
        case_by_key("gmres"),
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

    path = _reference_path(args, quadrature_order)
    with path.open("wb") as stream:
        pickle.dump(reference_data, stream, protocol=pickle.HIGHEST_PROTOCOL)
    print("wrote {}".format(path))
    return reference_data


def load_reference(args):
    path = _reference_path(args)
    if not path.exists():
        path = _legacy_reference_path(args)
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
            case_by_key("gs_si20"),
            tolerance,
            quadrature_order,
            args.pc_tolerance,
            args.linear_backend,
        )
        solver.solve()
        norm_phi, norm_fission = flux_norm(solver.state(), material, mesh, reference_data)
        rows.append((tolerance, norm_phi, norm_fission))
        print("{:12.5e} {:12.5e} {:12.5e}".format(tolerance, norm_phi, norm_fission))

    path = _tol_sweep_path(args, quadrature_order)
    np.savetxt(path, np.asarray(rows))
    print("wrote {}".format(path))


def run_convergence(args):
    quadrature_order = _quadrature(args)
    reference_data = load_reference(args)
    path = _convergence_path(args, quadrature_order)
    data = {
        "tolerance": float(args.tolerance),
        "quadrature": int(quadrature_order),
        "calibration_seed": args.calibration_seed,
        "tolerances": {},
        "sweeps": {},
    }
    if args.resume and path.exists():
        with path.open("rb") as stream:
            data = pickle.load(stream)
        data["calibration_seed"] = args.calibration_seed
        print("resuming from {}".format(path))
    tolerances = data["tolerances"]
    sweeps = data["sweeps"]

    for case in _selected_cases(args, CONVERGENCE_CASES):
        if args.resume and case.key in sweeps:
            print("{}: already calibrated; skipping".format(case.key))
            continue
        if args.calibration_seed == "paper":
            seed = _paper_residual_seed(args, case)
            if seed is not None:
                tolerance1 = float(seed)
                converged = False
                for iteration in range(50):
                    solver, material, mesh = get_solver(
                        case,
                        tolerance1,
                        quadrature_order,
                        args.pc_tolerance,
                        args.linear_backend,
                    )
                    solver.solve()
                    residual1, _ = flux_norm(
                        solver.state(),
                        material,
                        mesh,
                        reference_data,
                    )
                    print(
                        "{} paper-seed iter {:02d}: input tolerance {:12.5e}, error {:12.5e}".format(
                            case.key,
                            iteration,
                            tolerance1,
                            residual1,
                        )
                    )
                    sys.stdout.flush()
                    if residual1 < args.tolerance:
                        tolerances[case.code[0], case.code[1]] = tolerance1
                        sweeps[case.key] = solver.number_sweeps()
                        converged = True
                        break
                    tolerance1 *= 0.5
                if converged:
                    with path.open("wb") as stream:
                        pickle.dump(data, stream, protocol=pickle.HIGHEST_PROTOCOL)
                    print("checkpointed {}".format(path))
                    sys.stdout.flush()
                    continue
                print("{}: paper seed did not converge; falling back".format(case.key))
                sys.stdout.flush()
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
            sys.stdout.flush()
            if tolerance0 == tolerance1 or residual1 < args.tolerance:
                tolerances[case.code[0], case.code[1]] = tolerance1
                sweeps[case.key] = solver.number_sweeps()
                converged = True
                break
        if not converged:
            raise RuntimeError("could not calibrate {}".format(case.key))

        with path.open("wb") as stream:
            pickle.dump(data, stream, protocol=pickle.HIGHEST_PROTOCOL)
        print("checkpointed {}".format(path))
        sys.stdout.flush()

    with path.open("wb") as stream:
        pickle.dump(data, stream, protocol=pickle.HIGHEST_PROTOCOL)
    print("wrote {}".format(path))
    sys.stdout.flush()
    return data


def load_convergence(args):
    path_args = copy.copy(args)
    if args.calibration_source == "builtin":
        path_args.linear_backend = "builtin"
    path = _convergence_path(path_args)
    if not path.exists():
        path = _legacy_convergence_path(args)
    if not path.exists():
        raise FileNotFoundError(
            "missing {}; run `python3 run_all.py calibrate` first".format(path)
        )
    with path.open("rb") as stream:
        return pickle.load(stream)


def run_cases(args):
    quadrature_order = _quadrature(args)
    convergence = load_convergence(args)
    path = _case_data_path(args, quadrature_order)
    rows = []
    if args.resume and path.exists():
        with path.open("rb") as stream:
            rows = pickle.load(stream)
        print("resuming from {}".format(path))
    by_key = {row["case"]: row for row in rows}

    for case in _selected_cases(args, ARCHIVAL_CASES):
        if args.resume and case.key in by_key:
            print("{}: already run; skipping".format(case.key))
            continue
        if (case.code[0], case.code[1]) not in convergence["tolerances"]:
            raise RuntimeError(
                "missing calibrated tolerance for {}; run `python3 run_all.py "
                "calibrate --quadrature {} --tolerance {} --resume` first".format(
                    case.key,
                    quadrature_order,
                    _tol_label(args.tolerance),
                )
            )
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
        by_key[case.key] = rows[-1]
        print("{:26s} {:10.3f} s {:8d} sweeps".format(case.label, elapsed, sweeps))

        rows = _ordered_rows(by_key.values())
        with path.open("wb") as stream:
            pickle.dump(rows, stream, protocol=pickle.HIGHEST_PROTOCOL)
        print("checkpointed {}".format(path))

    rows = _ordered_rows(by_key.values())
    with path.open("wb") as stream:
        pickle.dump(rows, stream, protocol=pickle.HIGHEST_PROTOCOL)
    print("wrote {}".format(path))
    return rows


def load_case_data(args):
    quadrature_order = _quadrature(args)
    path = _case_data_path(args, quadrature_order)
    if not path.exists():
        raise FileNotFoundError(
            "missing {}; run `python3 run_all.py run` first".format(path)
        )
    with path.open("rb") as stream:
        return pickle.load(stream)


def print_table(args):
    rows = load_case_data(args)
    by_key = {row["case"]: row for row in rows}
    cases = _selected_cases(args, ARCHIVAL_CASES)
    missing = [case.key for case in cases if case.key not in by_key]
    if missing:
        raise RuntimeError(
            "missing case data for {}; rerun `python3 run_all.py run "
            "--quadrature {} --tolerance {} --resume` or select a smaller "
            "--cases list".format(
                ", ".join(missing),
                _quadrature(args),
                _tol_label(args.tolerance),
            )
        )

    out = []
    out.append("\\begin{table}[ht]")
    out.append("\\begin{center}")
    out.append("\\begin{tabular*}{0.7\\textwidth}{@{\\extracolsep{\\fill}} rccc}")
    out.append("\\toprule")
    out.append("                        & sweeps & sweeps/groups & time (s) \\\\")
    out.append("\\midrule")
    for case in cases:
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


def print_paper_table(args):
    paper = load_paper_results(args.paper_results)
    tol_key = _paper_tol_key(args.tolerance)
    tables = paper["c5g7_fixed_source"]["tables"]
    if tol_key not in tables:
        raise ValueError("paper table not available for tolerance {}".format(tol_key))

    print("Published C5G7 fixed-source table, tau = {}".format(tol_key))
    print("{:18s} {:>8s} {:>8s} {:>10s} {:>12s} {:>12s}".format(
        "case", "sweeps", "sweeps/g", "time(s)", "residual", "error"
    ))
    for row in tables[tol_key]:
        residual = "" if row["residual"] is None else "{:.3e}".format(row["residual"])
        print("{:18s} {:8d} {:8d} {:10.3f} {:>12s} {:12.3e}".format(
            row["label"],
            int(row["sweeps"]),
            int(row["sweeps_per_group"]),
            float(row["time_s"]),
            residual,
            float(row["error"]),
        ))


def compare_paper_table(args):
    paper = load_paper_results(args.paper_results)
    tol_key = _paper_tol_key(args.tolerance)
    paper_rows = {
        row["key"]: row
        for row in paper["c5g7_fixed_source"]["tables"].get(tol_key, [])
    }
    if not paper_rows:
        raise ValueError("paper table not available for tolerance {}".format(tol_key))

    try:
        rows = load_case_data(args)
    except FileNotFoundError as exc:
        raise SystemExit(str(exc))
    current_by_key = {row["case"]: row for row in rows}

    print(
        "Current {} C5G7 data vs published tau = {} table".format(
            args.linear_backend,
            tol_key,
        )
    )
    print("{:26s} {:>9s} {:>9s} {:>9s} {:>12s} {:>12s}".format(
        "current case", "now swp", "paper", "delta", "now time", "paper time"
    ))
    for case in _selected_cases(args, ARCHIVAL_CASES):
        paper_key = CURRENT_TO_PAPER_CASES.get(case.key)
        if paper_key is None or case.key not in current_by_key:
            continue
        current = current_by_key[case.key]
        published = paper_rows[paper_key]
        now_sweeps = int(current["sweeps"])
        paper_sweeps = int(published["sweeps"])
        delta = now_sweeps - paper_sweeps
        print("{:26s} {:9d} {:9d} {:9d} {:12.3f} {:12.3f}".format(
            case.label,
            now_sweeps,
            paper_sweeps,
            delta,
            float(current["time"]),
            float(published["time_s"]),
        ))
    if args.linear_backend == "petsc":
        print(
            "Note: paper timings are PETSc-era and hardware dependent; current "
            "PETSc timings should still be read primarily by sweep count."
        )
    else:
        print(
            "Note: paper timings used PETSc-era LU/ILU paths; built-in callow "
            "comparisons should be read primarily by sweep count."
        )


def smoke(args):
    quadrature_order = _quadrature(args, smoke=True)
    case = case_by_key("cmdsa")
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


def bwr_materials(paper):
    data = paper["bwr_spectral"]
    material = Material(len(data["materials"]), data["groups"], "paper-bwr")
    for entry in data["materials"]:
        m = int(entry["id"])
        for g, value in enumerate(entry["sigma_t"]):
            material.set_sigma_t(m, g, float(value))
        for g, value in enumerate(entry["sigma_f"]):
            material.set_sigma_f(m, g, float(value))
        for g, value in enumerate(entry["nu"]):
            material.set_nu(m, g, float(value))
        for g, value in enumerate(entry["chi"]):
            material.set_chi(m, g, float(value))
        for g_to, row in enumerate(entry["sigma_s"]):
            for g_from, value in enumerate(row):
                material.set_sigma_s(m, g_to, g_from, float(value))
    material.compute_diff_coef()
    material.finalize()
    return material


def bwr_coarse_map(paper, crop=9):
    interior = paper["bwr_spectral"]["interior_pin_map_bottom_to_top"]
    if crop < 1 or crop > len(interior) or crop % 2 == 0:
        raise ValueError("--bwr-crop must be an odd value from 1 through 9")
    start = (len(interior) - crop) // 2
    cropped = [row[start:start + crop] for row in interior[start:start + crop]]
    width = crop + 2
    moderator = [1] * width
    return [moderator] + [[1] + row + [1] for row in cropped] + [moderator]


def bwr_mesh(paper, crop=9, submesh=3):
    if submesh < 1:
        raise ValueError("--bwr-submesh must be positive")
    data = paper["bwr_spectral"]
    coarse_map = bwr_coarse_map(paper, crop)
    width = len(coarse_map)
    pitch = float(data["pitch_cm"])
    edges = [pitch * i for i in range(width + 1)]
    fine = [int(submesh)] * width
    mat_map = [int(value) for row in coarse_map for value in row]
    mesh = Mesh2D(fine, fine, edges, edges, mat_map)
    return mesh, coarse_map


def _bwr_solver_databases(args):
    return {
        "inner_solver_db": _linear_database(
            rtol=args.bwr_solver_tolerance,
            max_iters=max(20, args.bwr_max_iters),
            monitor_level=0,
            backend=args.linear_backend,
        ),
        "outer_solver_db": _linear_database(
            rtol=args.bwr_solver_tolerance,
            max_iters=args.bwr_max_iters,
            monitor_level=args.bwr_monitor_level,
            backend=args.linear_backend,
        ),
        "inner_pc_db": _preconditioner_database(
            "none",
            args.pc_tolerance,
            args.linear_backend,
        ),
        "outer_pc_db": _preconditioner_database(
            args.bwr_diffusion_pc,
            args.pc_tolerance,
            args.linear_backend,
        ),
    }


def bwr_overrides(args, kind):
    if kind in ("cmdsa", "cmdsa-s", "tcmdsa") and args.bwr_submesh < 2:
        raise ValueError(
            "CMDSA BWR operators need --bwr-submesh >= 2 so the coarse mesh "
            "is not identical to the fine mesh"
        )
    overrides = {
        "store_angular_flux": 0,
        "store_current": 0,
        "outer_krylov_group_cutoff": 0,
        "inner_pc_type": "none",
        "outer_pc_type": "none",
        "inner_pc_side": RIGHT,
        "outer_pc_side": RIGHT,
        "mgdsa_disable_fission": 0,
        "mgpc_disable_fission": int(args.bwr_disable_pc_fission),
        "mgpc_coarse_mesh_level": int(args.bwr_submesh),
        "mgpc_condensation_option": 1,
        "mgpc_cmdsa_use_smoothing": 0,
    }
    if kind == "mgdsa":
        overrides["outer_pc_type"] = "mgdsa"
    elif kind in ("cmdsa", "cmdsa-s", "tcmdsa"):
        overrides["outer_pc_type"] = "mgcmdsa" if kind in ("cmdsa", "cmdsa-s") else "mgtcdsa"
        overrides["mgpc_cmdsa_use_smoothing"] = int(kind in ("cmdsa-s", "tcmdsa"))
        if kind in ("cmdsa-s", "tcmdsa"):
            overrides["mgpc_cmdsa_smoothing_iters"] = 3
            overrides["mgpc_cmdsa_smoothing_relax"] = 0.7
        if kind == "tcmdsa":
            overrides["mgpc_tcdsa_number_coarse_corrections"] = 3
            overrides["mgpc_tcdsa_use_fine_correction"] = 1
            overrides["mgpc_number_azimuth_octant"] = 2
            overrides["mgpc_number_polar_octant"] = 2
    elif kind != "transport":
        raise ValueError("unsupported BWR operator kind {}".format(kind))

    return overrides


def bwr_input(args, kind):
    return fixed_source_input(
        groups=3,
        dimension=2,
        equation="dd",
        boundary="vacuum",
        inner=IterationSolver.source_iteration(
            max_iters=1,
            tolerance=1.0e-16,
            print_level=0,
        ),
        outer=IterationSolver.gmres(
            max_iters=args.bwr_max_iters,
            tolerance=args.bwr_solver_tolerance,
            print_level=args.bwr_monitor_level,
            print_interval=1,
        ),
        quadrature=Quadrature(
            type="u-dgl",
            polar=args.bwr_quadrature,
            azimuth=args.bwr_quadrature,
        ),
        solver_databases=_bwr_solver_databases(args),
        input_overrides=bwr_overrides(args, kind),
    )


def bwr_manager(args, paper, kind):
    material = bwr_materials(paper)
    mesh, coarse_map = bwr_mesh(paper, args.bwr_crop, args.bwr_submesh)
    problem = FixedSourceProblem(
        mesh=mesh,
        materials=material,
        boundary="vacuum",
        equation="dd",
        inner=IterationSolver.source_iteration(
            max_iters=1,
            tolerance=1.0e-16,
            print_level=0,
        ),
        outer=IterationSolver.gmres(
            max_iters=args.bwr_max_iters,
            tolerance=args.bwr_solver_tolerance,
            print_level=args.bwr_monitor_level,
            print_interval=1,
        ),
        quadrature=Quadrature(
            type="u-dgl",
            polar=args.bwr_quadrature,
            azimuth=args.bwr_quadrature,
        ),
        solver_databases=_bwr_solver_databases(args),
        input_overrides=bwr_overrides(args, kind),
        multiply=True,
        sources=[
            lambda manager: ConstantSource(
                3,
                manager.mesh(),
                1.0,
                manager.quadrature(),
            )
        ],
    )
    solver = problem.prepare_raw_solver()
    return solver, mesh, coarse_map


def _dense_gmres_iterations(A, tolerance=1.0e-8, max_iters=None):
    n = A.shape[0]
    if max_iters is None:
        max_iters = n
    b = np.ones(n)
    beta = np.linalg.norm(b)
    if beta == 0.0:
        return 0
    q = np.zeros((n, max_iters + 1))
    h = np.zeros((max_iters + 1, max_iters))
    q[:, 0] = b / beta
    e1 = np.zeros(max_iters + 1)
    e1[0] = beta
    for k in range(max_iters):
        v = A.dot(q[:, k])
        for j in range(k + 1):
            h[j, k] = np.dot(q[:, j], v)
            v -= h[j, k] * q[:, j]
        h[k + 1, k] = np.linalg.norm(v)
        if h[k + 1, k] != 0.0 and k + 1 < max_iters + 1:
            q[:, k + 1] = v / h[k + 1, k]
        rhs = e1[:k + 2]
        y, _, _, _ = np.linalg.lstsq(h[:k + 2, :k + 1], rhs, rcond=None)
        residual = np.linalg.norm(rhs - h[:k + 2, :k + 1].dot(y)) / beta
        if residual < tolerance:
            return k + 1
    return max_iters


def bwr_operator_summary(A, args):
    n, m = A.shape
    nnz = int(np.count_nonzero(np.abs(A) > args.drop_tolerance))
    summary = {
        "shape": [int(n), int(m)],
        "nnz": nnz,
        "density": float(nnz) / float(n * m),
    }
    if n == m:
        summary["frobenius_residual"] = float(np.linalg.norm(np.eye(n) - A, ord="fro"))
    if args.bwr_eigs:
        if n > args.bwr_max_diagnostics_size:
            print(
                "Skipping dense eigenvalue diagnostics for size {}; raise "
                "--bwr-max-diagnostics-size to force it.".format(n)
            )
        else:
            eigvals = np.linalg.eigvals(A)
            summary["std_lambda"] = float(np.std(eigvals))
            summary["spectral_radius"] = float(np.max(np.abs(eigvals)))
    if args.bwr_gmres:
        if n > args.bwr_max_diagnostics_size:
            print(
                "Skipping dense GMRES diagnostic for size {}; raise "
                "--bwr-max-diagnostics-size to force it.".format(n)
            )
        else:
            summary["gmres_iterations"] = int(
                _dense_gmres_iterations(A, 1.0e-8, min(n, args.bwr_max_iters))
            )
    return summary


def print_paper_bwr(args):
    paper = load_paper_results(args.paper_results)
    bwr = paper["bwr_spectral"]
    print("Published BWR spectral benchmark")
    print("operator size: {}".format(bwr["published_operator_size"]))
    print("nonzeros: {} ({:.3f} dense)".format(
        bwr["published_nonzeros"],
        bwr["published_density"],
    ))
    print("k_v ~= {}, k_inf ~= {}".format(bwr["vacuum_keff"], bwr["reflected_keff"]))
    print("{:10s} {:>12s} {:>14s} {:>10s}".format(
        "operator", "std(lambda)", "||I-AP^-1||F", "iters"
    ))
    for row in bwr["k1_spectrum_table"]:
        print("{:10s} {:12.4f} {:14.4f} {:10d}".format(
            row["key"],
            row["std_lambda"],
            row["frobenius_residual"],
            row["iterations"],
        ))


def run_bwr_operator(args):
    if args.bwr_paper_size:
        args.bwr_crop = 9
        args.bwr_submesh = 3
        args.bwr_quadrature = 6

    paper = load_paper_results(args.paper_results)
    solver, mesh, coarse_map = bwr_manager(args, paper, args.bwr_operator_kind)
    helper = getattr(_detran, "mg_transport_operator_dense_2d", None)
    size_helper = getattr(_detran, "mg_transport_operator_size_2d", None)
    if helper is None or size_helper is None:
        raise RuntimeError(
            "pydetran does not expose mg_transport_operator_dense_2d; rebuild "
            "and reinstall pydetran from the current branch"
        )

    shape = list(size_helper(solver))
    print(
        "BWR {} operator: shape {}x{}, cells {}, crop {}, submesh {}, q {}".format(
            args.bwr_operator_kind,
            shape[0],
            shape[1],
            mesh.number_cells(),
            args.bwr_crop,
            args.bwr_submesh,
            args.bwr_quadrature,
        )
    )
    if args.bwr_size_only:
        return {"shape": shape, "coarse_map": coarse_map}

    right_preconditioned = args.bwr_operator_kind != "transport"
    start = time.time()
    flat = helper(
        solver,
        args.keff,
        right_preconditioned,
        args.drop_tolerance,
    )
    elapsed = time.time() - start
    A = np.asarray(flat, dtype=float).reshape((shape[1], shape[0])).T
    nonfinite = int(np.size(A) - np.count_nonzero(np.isfinite(A)))
    if nonfinite:
        raise RuntimeError(
            "{} BWR operator contains {} non-finite entries; try a different "
            "--keff, --bwr-submesh, or --bwr-diffusion-pc setting".format(
                args.bwr_operator_kind,
                nonfinite,
            )
        )
    summary = bwr_operator_summary(A, args)
    summary.update({
        "operator_kind": args.bwr_operator_kind,
        "keff": args.keff,
        "elapsed_s": elapsed,
        "crop": args.bwr_crop,
        "submesh": args.bwr_submesh,
        "quadrature": args.bwr_quadrature,
        "pc_fission": not bool(args.bwr_disable_pc_fission),
    })

    output = Path(args.output) if args.output else _artifact(
        args.data_dir,
        "bwr_operator_{}_{}_{}_q{}_k{}.npz".format(
            args.bwr_operator_kind,
            args.bwr_crop,
            args.bwr_submesh,
            args.bwr_quadrature,
            _tol_label(args.keff),
        ),
    )
    np.savez_compressed(
        output,
        operator=A,
        summary=json.dumps(summary, sort_keys=True),
        coarse_map=np.asarray(coarse_map, dtype=int),
    )
    print("assembled in {:.3f} s".format(elapsed))
    print(json.dumps(summary, indent=2, sort_keys=True))
    print("wrote {}".format(output))
    return summary


def run_bwr_operator_suite(args, *, paper_size=False, dense=False, diagnostics=True):
    summaries = {}
    kinds = BWR_OPERATOR_KINDS if dense or not paper_size else ("transport",)
    for kind in kinds:
        local = copy.copy(args)
        local.bwr_operator_kind = kind
        local.bwr_paper_size = bool(paper_size)
        local.bwr_size_only = bool(paper_size and not dense)
        local.bwr_eigs = bool(diagnostics and (dense or not paper_size))
        local.bwr_gmres = bool(diagnostics and (dense or not paper_size))
        summaries[kind] = run_bwr_operator(local)
    return summaries


def _paper_manifest_path(args):
    return _artifact(
        args.data_dir,
        "paper_reproduction_{}.json".format(args.paper_mode),
    )


def _write_paper_manifest(args, manifest):
    path = _paper_manifest_path(args)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print("wrote {}".format(path))
    return path


def print_paper_plan(args):
    print("Paper reproduction coverage")
    print("  C5G7 tables: reference, tolerance calibration, runs, table, compare")
    print("  C5G7 norms: tolerance-sweep data for flux vs pin-power error")
    print("  BWR spectra: reduced dense operators by default")
    print("  BWR full-size: size check by default; dense assembly is opt-in")
    print("")
    print("Quick smoke:")
    print("  python3 run_all.py reproduce-paper --paper-mode smoke")
    print("")
    print("One current tolerance:")
    print("  python3 run_all.py reproduce-paper --paper-mode current --quadrature 8 --tolerance 1e-6")
    print("")
    print("Full C5G7 paper tables:")
    print("  python3 run_all.py reproduce-paper --paper-mode full --quadrature 8")
    print("")
    print("Full BWR dense operators, intentionally expensive:")
    print("  python3 run_all.py reproduce-paper --paper-mode full --include-full-bwr-dense")


def reproduce_paper(args):
    paper = load_paper_results(args.paper_results)
    manifest = {
        "paper": paper["paper"],
        "mode": args.paper_mode,
        "linear_backend": args.linear_backend,
        "quadrature": _quadrature(args, smoke=args.paper_mode == "smoke"),
        "c5g7": {},
        "bwr": {},
        "limitations": [
            "TC^2-CMDSA+S is imported from the paper table but is not yet "
            "implemented as a current run case.",
            "Default built-in LU rows use full-fill ILUT as the closest "
            "available callow equivalent; PETSc mode uses PETSc LU.",
            "Paper timings are hardware and PETSc-era dependent; current "
            "comparisons should be read primarily by sweep count.",
        ],
    }

    if args.paper_mode == "plan":
        print_paper_plan(args)
        manifest["status"] = "plan-only"
        _write_paper_manifest(args, manifest)
        return manifest

    if args.paper_mode == "smoke":
        smoke(args)
        manifest["c5g7"]["smoke"] = {
            "case": "cmdsa",
            "quadrature": _quadrature(args, smoke=True),
        }
        manifest["bwr"]["reduced"] = run_bwr_operator_suite(
            args,
            paper_size=False,
            dense=True,
            diagnostics=True,
        )
        size_args = copy.copy(args)
        size_args.bwr_paper_size = True
        size_args.bwr_size_only = True
        size_args.bwr_operator_kind = "transport"
        manifest["bwr"]["paper_size"] = run_bwr_operator(size_args)
        _write_paper_manifest(args, manifest)
        return manifest

    if args.paper_mode == "current":
        run_reference(args)
        if not args.skip_tolerance_sweep:
            flux_norm_vs_pin_power(args)
        run_convergence(args)
        rows = run_cases(args)
        print_table(args)
        compare_paper_table(args)
        manifest["c5g7"][_paper_tol_key(args.tolerance)] = rows
        manifest["bwr"]["reduced"] = run_bwr_operator_suite(
            args,
            paper_size=False,
            dense=True,
            diagnostics=True,
        )
        _write_paper_manifest(args, manifest)
        return manifest

    if args.paper_mode == "full":
        run_reference(args)
        if not args.skip_tolerance_sweep:
            flux_norm_vs_pin_power(args)
        for tolerance in args.paper_tolerances:
            local = copy.copy(args)
            local.tolerance = float(tolerance)
            run_convergence(local)
            rows = run_cases(local)
            print_table(local)
            compare_paper_table(local)
            manifest["c5g7"][_paper_tol_key(tolerance)] = rows
        manifest["bwr"]["paper_size"] = run_bwr_operator_suite(
            args,
            paper_size=True,
            dense=bool(args.include_full_bwr_dense),
            diagnostics=bool(args.include_full_bwr_dense),
        )
        _write_paper_manifest(args, manifest)
        return manifest

    raise ValueError("unsupported paper mode {}".format(args.paper_mode))


def _relative_to_script(path):
    path = Path(path)
    try:
        return str(path.relative_to(SCRIPT_DIR))
    except ValueError:
        return str(path)


def _json_case_rows(rows):
    return [
        {
            "case": str(row["case"]),
            "backend": str(row.get("backend", "")),
            "sweeps": int(row["sweeps"]),
            "time_s": float(row["time"]),
        }
        for row in rows
    ]


def _json_convergence(data):
    if data is None:
        return None
    tolerances = []
    for code, tolerance in sorted(data.get("tolerances", {}).items()):
        tolerances.append({
            "solver_code": [int(value) for value in code],
            "tolerance": float(tolerance),
        })
    return {
        "tolerance": float(data["tolerance"]),
        "quadrature": int(data["quadrature"]),
        "calibration_seed": data.get("calibration_seed", "default"),
        "solver_tolerances": tolerances,
        "calibration_sweeps": {
            str(key): int(value)
            for key, value in sorted(data.get("sweeps", {}).items())
        },
    }


def _load_pickle_if_present(path):
    path = Path(path)
    if not path.exists():
        return None
    with path.open("rb") as stream:
        return pickle.load(stream)


def _archive_c5g7_entry(*, label, backend, tolerance, quadrature, rows_path, convergence_path):
    rows = _load_pickle_if_present(rows_path)
    if rows is None:
        return None
    convergence = _load_pickle_if_present(convergence_path)
    sources = [_relative_to_script(rows_path)]
    if Path(convergence_path).exists():
        sources.append(_relative_to_script(convergence_path))
    return {
        "label": label,
        "backend": backend,
        "tolerance": float(tolerance),
        "tolerance_key": _paper_tol_key(tolerance),
        "quadrature": int(quadrature),
        "source_files": sources,
        "rows": _json_case_rows(rows),
        "convergence": _json_convergence(convergence),
    }


def _write_updated_results_readme(output_dir):
    text = """# Updated Local Results

This directory is a compact archive for the local CMDSA reproduction results.
It stores only plotting inputs and the generated PDF report; dense BWR operator
matrices and raw log files are intentionally left as ignored generated data.

Files:

- `c5g7_local_runs.json`: local C5G7 timing and sweep rows.
- `bwr_spectral_diagnostics_with_fission.json`: BWR current-vs-paper metrics.
- `bwr_spectral_eigenvalues_with_fission.npz`: BWR eigenvalue arrays used for
  spectral plots.
- `updated_results.pdf`: standalone report generated from the archived data.

Refresh from local generated artifacts:

```bash
python3 ../run_all.py archive-results
python3 ../run_all.py updated-results
```
"""
    path = Path(output_dir) / "README.md"
    with path.open("w", encoding="utf-8") as stream:
        stream.write(text)
    return path


def archive_updated_results(args):
    output_dir = Path(args.updated_results_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paper = load_paper_results(args.paper_results)

    c5g7_specs = [
        {
            "label": "builtin q8 tau=1e-4",
            "backend": "builtin",
            "tolerance": 1.0e-4,
            "quadrature": 8,
            "rows_path": SCRIPT_DIR / "Xdata_0.0001_8.p",
            "convergence_path": SCRIPT_DIR / "convergence_data_0.0001_q8.p",
        },
        {
            "label": "petsc 1thread q8 tau=1e-4",
            "backend": "petsc",
            "tolerance": 1.0e-4,
            "quadrature": 8,
            "rows_path": SCRIPT_DIR / "thread_pinning_runs" / "petsc_1thread" / "Xdata_0.0001_8_petsc.p",
            "convergence_path": SCRIPT_DIR / "thread_pinning_runs" / "petsc_1thread" / "convergence_data_0.0001_q8.p",
        },
        {
            "label": "petsc 1thread q8 tau=1e-6",
            "backend": "petsc",
            "tolerance": 1.0e-6,
            "quadrature": 8,
            "rows_path": SCRIPT_DIR / "thread_pinning_runs" / "petsc_1thread" / "Xdata_1e-06_8_petsc.p",
            "convergence_path": SCRIPT_DIR / "thread_pinning_runs" / "petsc_1thread" / "convergence_data_1e-06_q8_petsc.p",
        },
        {
            "label": "petsc 1thread q8 tau=1e-8",
            "backend": "petsc",
            "tolerance": 1.0e-8,
            "quadrature": 8,
            "rows_path": SCRIPT_DIR / "thread_pinning_runs" / "petsc_1thread" / "Xdata_1e-08_8_petsc.p",
            "convergence_path": SCRIPT_DIR / "thread_pinning_runs" / "petsc_1thread" / "convergence_data_1e-08_q8_petsc.p",
        },
    ]
    runs = []
    for spec in c5g7_specs:
        entry = _archive_c5g7_entry(**spec)
        if entry is not None:
            runs.append(entry)

    c5g7_archive = {
        "description": "Compact local C5G7 timing and sweep data for updated results plotting.",
        "paper": paper["paper"],
        "runs": runs,
    }
    c5g7_path = output_dir / UPDATED_RESULTS_C5G7
    with c5g7_path.open("w", encoding="utf-8") as stream:
        json.dump(c5g7_archive, stream, indent=2, sort_keys=True)
        stream.write("\n")

    spectral_dir = Path(args.spectral_results_dir)
    diagnostics_path = spectral_dir / "spectral_diagnostics_with_fission.json"
    if not diagnostics_path.exists():
        raise FileNotFoundError(
            "missing {}; regenerate BWR spectral diagnostics first".format(
                diagnostics_path,
            )
        )
    with diagnostics_path.open("r", encoding="utf-8") as stream:
        bwr_diagnostics = json.load(stream)
    bwr_diagnostics["_archive_notes"] = [
        "Current A/DSA/CMDSA/CMDSA-S entries use corrected BWR fission data.",
        "TC-CMDSA was not regenerated in corrected-fission mode for this archive.",
        "Dense operator matrices are excluded; eigenvalues and diagnostics are enough to replot.",
    ]
    bwr_path = output_dir / UPDATED_RESULTS_BWR_DIAGNOSTICS
    with bwr_path.open("w", encoding="utf-8") as stream:
        json.dump(bwr_diagnostics, stream, indent=2, sort_keys=True)
        stream.write("\n")

    eig_arrays = {}
    for label, (archive_key, filename) in BWR_EIGENVALUE_SOURCES.items():
        source = spectral_dir / filename
        if not source.exists():
            raise FileNotFoundError("missing {}".format(source))
        with np.load(source) as data:
            eig_arrays[archive_key] = data["eigvals"]
    eig_path = output_dir / UPDATED_RESULTS_BWR_EIGS
    np.savez_compressed(eig_path, **eig_arrays)

    readme_path = _write_updated_results_readme(output_dir)
    print("wrote {}".format(c5g7_path))
    print("wrote {}".format(bwr_path))
    print("wrote {}".format(eig_path))
    print("wrote {}".format(readme_path))
    return {
        "c5g7": str(c5g7_path),
        "bwr": str(bwr_path),
        "eigs": str(eig_path),
        "readme": str(readme_path),
    }


def _case_label(case_key):
    for case in ARCHIVAL_CASES:
        if case.key == case_key:
            return case.label
    return case_key


def _paper_rows_by_key(paper, tolerance_key):
    return {
        row["key"]: row
        for row in paper["c5g7_fixed_source"]["tables"].get(tolerance_key, [])
    }


def _report_case_rows(run, paper):
    paper_rows = _paper_rows_by_key(paper, run["tolerance_key"])
    rows = []
    for row in run["rows"]:
        paper_key = CURRENT_TO_PAPER_CASES.get(row["case"])
        if paper_key not in paper_rows:
            continue
        published = paper_rows[paper_key]
        rows.append({
            "label": _case_label(row["case"]),
            "current_sweeps": int(row["sweeps"]),
            "paper_sweeps": int(published["sweeps"]),
            "current_time": float(row["time_s"]),
            "paper_time": float(published["time_s"]),
        })
    return rows


def _short_label(label):
    return (
        label.replace("GMRES+", "")
        .replace("GS, ", "GS ")
        .replace("TC-CMDSA", "TC")
        .replace("CMDSA", "CM")
    )


def _load_updated_results(args):
    result_dir = Path(args.updated_results_dir)
    with (result_dir / UPDATED_RESULTS_C5G7).open("r", encoding="utf-8") as stream:
        c5g7 = json.load(stream)
    with (result_dir / UPDATED_RESULTS_BWR_DIAGNOSTICS).open("r", encoding="utf-8") as stream:
        bwr = json.load(stream)
    eigs = np.load(result_dir / UPDATED_RESULTS_BWR_EIGS)
    return result_dir, c5g7, bwr, eigs


def write_updated_results_report(args):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required to write the updated-results PDF"
        ) from exc

    result_dir, c5g7, bwr, eigs = _load_updated_results(args)
    paper = load_paper_results(args.paper_results)
    output = Path(args.updated_results_pdf) if args.updated_results_pdf else (
        result_dir / UPDATED_RESULTS_REPORT
    )

    with PdfPages(output) as pdf:
        fig = plt.figure(figsize=(11.0, 8.5))
        fig.text(0.08, 0.86, "CMDSA Python Example: Updated Local Results", fontsize=20)
        paper_meta = paper["paper"]
        fig.text(
            0.08,
            0.78,
            "{}\n{} and {}, JCP {} ({}), pp. {}\n{}".format(
                paper_meta["title"],
                paper_meta["authors"][0],
                paper_meta["authors"][1],
                paper_meta["volume"],
                paper_meta["year"],
                paper_meta["pages"],
                paper_meta["url"],
            ),
            fontsize=11,
            linespacing=1.35,
        )
        notes = [
            "Archive source: {}".format(result_dir),
            "C5G7 rows are local q8 timing/sweep runs compared with the recovered paper tables.",
            "BWR spectral rows use corrected fission material data for A/DSA/CMDSA/CMDSA-S.",
            "Corrected-fission TC-CMDSA was not regenerated in this pass.",
            "Transport nnz uses the archive drop tolerance; tiny entries near 1e-9 explain most nnz drift.",
        ]
        fig.text(0.08, 0.52, "\n".join("- " + note for note in notes), fontsize=11)
        pdf.savefig(fig)
        plt.close(fig)

        runs = c5g7["runs"]
        cols = 2
        rows_count = int(np.ceil(float(len(runs)) / cols))
        fig, axes = plt.subplots(rows_count, cols, figsize=(11.0, 8.5), squeeze=False)
        for ax, run in zip(axes.ravel(), runs):
            rows = _report_case_rows(run, paper)
            labels = [_short_label(row["label"]) for row in rows]
            y = np.arange(len(rows))
            ax.barh(y - 0.18, [row["paper_sweeps"] for row in rows], height=0.35, label="paper")
            ax.barh(y + 0.18, [row["current_sweeps"] for row in rows], height=0.35, label="local")
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7)
            ax.invert_yaxis()
            ax.set_xlabel("sweeps")
            ax.set_title(run["label"], fontsize=10)
            ax.grid(True, axis="x", alpha=0.25)
        for ax in axes.ravel()[len(runs):]:
            ax.axis("off")
        axes[0, 0].legend(loc="lower right", fontsize=8)
        fig.suptitle("C5G7 Sweep Counts", fontsize=15)
        fig.tight_layout(rect=(0.02, 0.02, 1.0, 0.95))
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(rows_count, cols, figsize=(11.0, 8.5), squeeze=False)
        for ax, run in zip(axes.ravel(), runs):
            rows = _report_case_rows(run, paper)
            labels = [_short_label(row["label"]) for row in rows]
            ratios = [
                row["current_time"] / row["paper_time"]
                if row["paper_time"] > 0.0 else np.nan
                for row in rows
            ]
            y = np.arange(len(rows))
            ax.barh(y, ratios, color="#4c78a8")
            ax.axvline(1.0, color="black", linewidth=0.8)
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7)
            ax.invert_yaxis()
            ax.set_xlabel("local time / paper time")
            ax.set_title(run["label"], fontsize=10)
            ax.grid(True, axis="x", alpha=0.25)
        for ax in axes.ravel()[len(runs):]:
            ax.axis("off")
        fig.suptitle("C5G7 Timing Ratios", fontsize=15)
        fig.tight_layout(rect=(0.02, 0.02, 1.0, 0.95))
        pdf.savefig(fig)
        plt.close(fig)

        labels = []
        paper_std = []
        current_std = []
        paper_fro = []
        current_fro = []
        for label in BWR_REPORT_ORDER:
            entry = bwr.get(label, {})
            if "paper" not in entry:
                continue
            labels.append(label)
            paper_std.append(float(entry["paper"]["std_lambda"]))
            paper_fro.append(float(entry["paper"]["frobenius_residual"]))
            current = entry.get("current")
            current_std.append(np.nan if current is None else float(current["std_lambda"]))
            current_fro.append(np.nan if current is None else float(current["frobenius_residual"]))
        x = np.arange(len(labels))
        width = 0.36
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.2))
        axes[0].bar(x - width / 2, paper_std, width, label="paper")
        axes[0].bar(x + width / 2, current_std, width, label="local")
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(labels, rotation=25, ha="right")
        axes[0].set_ylabel("std(lambda)")
        axes[0].grid(True, axis="y", alpha=0.25)
        axes[0].legend()
        axes[1].bar(x - width / 2, paper_fro, width, label="paper")
        axes[1].bar(x + width / 2, current_fro, width, label="local")
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(labels, rotation=25, ha="right")
        axes[1].set_ylabel("||I - A P^-1||_F")
        axes[1].grid(True, axis="y", alpha=0.25)
        fig.suptitle("BWR Spectral Metrics, k = 1", fontsize=15)
        fig.tight_layout(rect=(0.02, 0.02, 1.0, 0.92))
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.5), squeeze=False)
        for ax, label in zip(axes.ravel(), ("A", "DSA", "CMDSA", "CMDSA-S")):
            archive_key = BWR_EIGENVALUE_SOURCES[label][0]
            values = eigs[archive_key]
            ax.scatter(values.real, values.imag, s=3, alpha=0.35)
            ax.axhline(0.0, color="black", linewidth=0.5)
            ax.axvline(0.0, color="black", linewidth=0.5)
            current = bwr[label]["current"]
            ax.set_title(
                "{}: std={:.4f}, rho={:.4f}".format(
                    label,
                    float(current["std_lambda"]),
                    float(current["spectral_radius"]),
                ),
                fontsize=10,
            )
            ax.set_xlabel("real(lambda)")
            ax.set_ylabel("imag(lambda)")
            ax.grid(True, alpha=0.2)
        fig.suptitle("BWR Eigenvalue Spectra From Archived Local Run", fontsize=15)
        fig.tight_layout(rect=(0.02, 0.02, 1.0, 0.95))
        pdf.savefig(fig)
        plt.close(fig)

    eigs.close()
    print("wrote {}".format(output))
    return output


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        nargs="?",
        default="smoke",
        choices=(
            "smoke",
            "reference",
            "tolerance-sweep",
            "calibrate",
            "run",
            "table",
            "all",
            "paper-table",
            "paper-compare",
            "paper-bwr",
            "bwr-operator",
            "reproduce-paper",
            "archive-results",
            "updated-results",
        ),
    )
    parser.add_argument("--quadrature", type=int)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument("--pc-tolerance", type=float, default=0.01)
    parser.add_argument(
        "--cases",
        nargs="+",
        help="case keys to calibrate, run, print, or compare; defaults to all",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="reuse completed calibration or timing rows from an existing artifact",
    )
    parser.add_argument("--paper-results", default=str(PAPER_RESULTS))
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
    parser.add_argument(
        "--calibration-seed",
        choices=("default", "paper"),
        default="default",
        help="initial guesses for calibrate; paper verifies published residuals before accepting them",
    )
    parser.add_argument("--data-dir", default=str(SCRIPT_DIR))
    parser.add_argument("--output")
    parser.add_argument("--updated-results-dir", default=str(UPDATED_RESULTS_DIR))
    parser.add_argument("--updated-results-pdf")
    parser.add_argument("--spectral-results-dir", default=str(SCRIPT_DIR / "spectral_plots"))
    parser.add_argument(
        "--paper-mode",
        choices=("plan", "smoke", "current", "full"),
        default="smoke",
        help="scope for reproduce-paper",
    )
    parser.add_argument(
        "--paper-tolerances",
        nargs="+",
        type=float,
        default=[1.0e-4, 1.0e-6, 1.0e-8],
        help="true tolerances for reproduce-paper --paper-mode full",
    )
    parser.add_argument(
        "--skip-tolerance-sweep",
        action="store_true",
        help="skip C5G7 flux-vs-pin-power tolerance-sweep data",
    )
    parser.add_argument(
        "--include-full-bwr-dense",
        action="store_true",
        help="assemble full paper-size dense BWR operators in full mode",
    )
    parser.add_argument("--keff", type=float, default=1.0)
    parser.add_argument("--drop-tolerance", type=float, default=1.0e-14)
    parser.add_argument("--bwr-crop", type=int, default=1)
    parser.add_argument("--bwr-submesh", type=int, default=3)
    parser.add_argument("--bwr-quadrature", type=int, default=1)
    parser.add_argument("--bwr-paper-size", action="store_true")
    parser.add_argument("--bwr-size-only", action="store_true")
    parser.add_argument(
        "--bwr-operator-kind",
        choices=BWR_OPERATOR_KINDS,
        default="transport",
    )
    parser.add_argument("--bwr-eigs", action="store_true")
    parser.add_argument("--bwr-gmres", action="store_true")
    parser.add_argument("--bwr-max-diagnostics-size", type=int, default=300)
    parser.add_argument("--bwr-max-iters", type=int, default=5000)
    parser.add_argument("--bwr-monitor-level", type=int, default=0)
    parser.add_argument("--bwr-solver-tolerance", type=float, default=1.0e-10)
    parser.add_argument(
        "--bwr-diffusion-pc",
        choices=("none", "ilu0", "ilut"),
        default="none",
        help="preconditioner used inside BWR diffusion acceleration solves",
    )
    parser.add_argument(
        "--bwr-disable-pc-fission",
        action="store_true",
        help="omit fission from BWR CMDSA-family preconditioner diagnostics",
    )
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
    elif args.action == "paper-table":
        print_paper_table(args)
    elif args.action == "paper-compare":
        compare_paper_table(args)
    elif args.action == "paper-bwr":
        print_paper_bwr(args)
    elif args.action == "bwr-operator":
        run_bwr_operator(args)
    elif args.action == "reproduce-paper":
        reproduce_paper(args)
    elif args.action == "archive-results":
        archive_updated_results(args)
    elif args.action == "updated-results":
        write_updated_results_report(args)
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
