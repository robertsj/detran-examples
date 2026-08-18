"""Run centered-mesh MOC pin-cell refinement sweeps.

The default sweep is intentionally pin-cell only.  It doubles the centered
Cartesian pin mesh from 4x4 upward, uses uniform azimuthal angles with
half-space Gauss-Legendre polar quadrature (`u-gl`), and stops before launching
the next run when the next wall time is predicted to exceed the requested
limit.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Iterable

os.environ.setdefault("OMP_NUM_THREADS", "16")

import numpy as np
import pandas as pd

try:
  from .moc_driver import (
    MOC_PROBLEM_OPTION_NAMES,
    build_moc_problem,
    get_moc_boundary,
    get_moc_case,
    material_for_mesh_recipe,
    summarize_moc_result,
  )
except ImportError:
  from moc_driver import (
    MOC_PROBLEM_OPTION_NAMES,
    build_moc_problem,
    get_moc_boundary,
    get_moc_case,
    material_for_mesh_recipe,
    summarize_moc_result,
  )


@dataclass(frozen=True)
class SweepRun:
  run_id: str
  case: str
  index: int
  mesh_n: int
  cells: int
  azimuth: int
  polar: int
  quadrature: str
  spacing: float
  spacing_mode: str
  density_correction: bool
  equation: str
  elapsed_s: float | None
  predicted_next_s: float | None
  stop_reason: str | None
  status: str
  keff: float | None = None
  finite: bool | None = None
  sweeps: int | None = None
  reference_keff: float | None = None
  delta_keff: float | None = None
  error: str | None = None
  result_path: str | None = None


def _parse_int_list(spec: str) -> list[int]:
  values = []
  for part in str(spec).split(","):
    part = part.strip()
    if part:
      values.append(int(part))
  if not values:
    raise argparse.ArgumentTypeError("list must contain at least one integer")
  return values


def _mesh_values(start: int, growth: int, max_mesh: int) -> Iterable[int]:
  current = int(start)
  while current <= int(max_mesh):
    yield current
    current *= int(growth)


def _cell_spacing(mesh) -> float:
  dx = [float(mesh.dx(i)) for i in range(int(mesh.number_cells_x()))]
  dy = [float(mesh.dy(j)) for j in range(int(mesh.number_cells_y()))]
  return min(dx + dy)


def _selected_spacing(mesh, args) -> float:
  if args.spacing_mode == "uniform":
    return float(args.spacing)
  if args.spacing_mode == "cell":
    return float(args.spacing_scale) * _cell_spacing(mesh)
  raise ValueError("spacing_mode must be uniform or cell")


def _mesh_materials(mesh) -> np.ndarray:
  if not mesh.mesh_map_exists("MATERIAL"):
    raise RuntimeError("mesh does not expose a MATERIAL map")
  return np.asarray(mesh.mesh_map("MATERIAL"), dtype=np.int64)


def _cell_indices(mesh, cell: int) -> tuple[int, int]:
  nx = int(mesh.number_cells_x())
  return int(cell % nx), int(cell // nx)


def _flux_tables(run_id: str, mesh, material, result) -> tuple[pd.DataFrame, pd.DataFrame]:
  number_cells = int(mesh.number_cells())
  number_groups = int(material.number_groups())
  materials = _mesh_materials(mesh)
  volumes = np.asarray([float(mesh.volume(c)) for c in range(number_cells)], dtype=float)
  cell_rows = []
  for group in range(number_groups):
    flux = np.asarray(result.flux(group), dtype=float)
    for cell in range(number_cells):
      i, j = _cell_indices(mesh, cell)
      cell_rows.append({
        "run_id": run_id,
        "cell": cell,
        "i": i,
        "j": j,
        "material": int(materials[cell]),
        "group": group,
        "volume": float(volumes[cell]),
        "flux": float(flux[cell]),
      })
  cell_df = pd.DataFrame(cell_rows)

  material_rows = []
  for group in range(number_groups):
    flux = np.asarray(result.flux(group), dtype=float)
    for mat in sorted(int(value) for value in np.unique(materials)):
      mask = materials == mat
      volume = float(np.sum(volumes[mask]))
      mean_flux = float(np.sum(flux[mask] * volumes[mask]) / volume)
      material_rows.append({
        "run_id": run_id,
        "material": mat,
        "group": group,
        "volume": volume,
        "mean_flux": mean_flux,
      })
  return cell_df, pd.DataFrame(material_rows)


def _material_area_rows(run_id: str, recipe: dict) -> list[dict]:
  rows = []
  for entry in recipe.get("material_area_corrections", []):
    row = {"run_id": run_id}
    row.update(entry)
    rows.append(row)
  return rows


def _git_commit() -> str | None:
  try:
    return subprocess.check_output(
      ["git", "rev-parse", "HEAD"],
      text=True,
      stderr=subprocess.DEVNULL,
    ).strip()
  except Exception:
    return None


def _run_case(args, mesh_n: int, polar: int, azimuth: int, output_dir: Path) -> tuple[
  SweepRun,
  pd.DataFrame,
  pd.DataFrame,
  list[dict],
  float | None,
]:
  correction_label = "corr" if args.density_correction else "uncorr"
  run_id = (
    f"pin{args.index}_n{mesh_n}_p{polar}_a{azimuth}"
    f"_s{args.spacing_mode}_{correction_label}"
  )
  run_dir = output_dir / "runs"
  run_dir.mkdir(parents=True, exist_ok=True)
  result_path = run_dir / f"{run_id}.json"

  options = {
    "index": args.index,
    "mesh_variant": "centered",
    "pin_mesh": mesh_n,
    "density_correction": args.density_correction,
    "polar": polar,
    "azimuth": azimuth,
    "quadrature": args.quadrature,
    "spacing": None,
    "equation": args.equation,
    "inner_tolerance": args.inner_tolerance,
    "outer_tolerance": args.outer_tolerance,
    "eigen_tolerance": args.eigen_tolerance,
    "inner_iters": args.inner_iters,
    "outer_iters": args.outer_iters,
    "eigen_iters": args.eigen_iters,
    "monitor_level": args.monitor_level,
    "print_interval": args.print_interval,
    "plot_dir": None,
    "boundary": get_moc_boundary(args.case),
  }

  mesh, mesh_recipe, reference = get_moc_case(
    args.case,
    index=args.index,
    pin_mesh=mesh_n,
    mesh_variant="centered",
  )
  options["spacing"] = _selected_spacing(mesh, args)
  material, material_model = material_for_mesh_recipe(
    mesh_recipe,
    apply_density_correction=args.density_correction,
  )
  problem_options = {name: options[name] for name in MOC_PROBLEM_OPTION_NAMES}
  problem_options["boundary"] = options["boundary"]
  problem = build_moc_problem(mesh, material, **problem_options)

  start = time.time()
  try:
    result = problem.solve()
    elapsed = time.time() - start
    summary = summarize_moc_result(
      args.case,
      args.index,
      mesh,
      result,
      reference,
      options,
      elapsed,
      mesh_recipe,
      material_model,
    )
    summary.update({
      "run_id": run_id,
      "status": "completed",
      "spacing_mode": args.spacing_mode,
      "spacing_scale": args.spacing_scale,
      "threads": args.threads,
      "solver_model": "flattened_power_iteration",
    })
    cell_df, material_df = _flux_tables(run_id, mesh, material, result)
    area_rows = _material_area_rows(run_id, mesh_recipe)
    result_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    run = SweepRun(
      run_id=run_id,
      case=args.case,
      index=args.index,
      mesh_n=mesh_n,
      cells=int(mesh.number_cells()),
      azimuth=azimuth,
      polar=polar,
      quadrature=args.quadrature,
      spacing=float(options["spacing"]),
      spacing_mode=args.spacing_mode,
      density_correction=bool(args.density_correction),
      equation=args.equation,
      elapsed_s=elapsed,
      predicted_next_s=None,
      stop_reason=None,
      status="completed",
      keff=float(summary["keff"]),
      finite=bool(summary["finite"]),
      sweeps=int(summary.get("sweeps", 0)),
      reference_keff=float(summary.get("reference_keff")),
      delta_keff=float(summary.get("delta_keff")),
      result_path=str(result_path),
    )
    return run, cell_df, material_df, area_rows, elapsed
  except Exception as exc:
    elapsed = time.time() - start
    failure = {
      "run_id": run_id,
      "status": "failed",
      "error": repr(exc),
      "case": args.case,
      "index": args.index,
      "mesh_n": mesh_n,
      "polar": polar,
      "azimuth": azimuth,
      "quadrature": args.quadrature,
      "spacing": float(options["spacing"]),
      "mesh_recipe": mesh_recipe,
      "material_model": material_model,
      "elapsed_s": elapsed,
    }
    result_path.write_text(json.dumps(failure, indent=2, sort_keys=True) + "\n")
    run = SweepRun(
      run_id=run_id,
      case=args.case,
      index=args.index,
      mesh_n=mesh_n,
      cells=int(mesh.number_cells()),
      azimuth=azimuth,
      polar=polar,
      quadrature=args.quadrature,
      spacing=float(options["spacing"]),
      spacing_mode=args.spacing_mode,
      density_correction=bool(args.density_correction),
      equation=args.equation,
      elapsed_s=elapsed,
      predicted_next_s=None,
      stop_reason="failed",
      status="failed",
      error=repr(exc),
      result_path=str(result_path),
    )
    return run, pd.DataFrame(), pd.DataFrame(), _material_area_rows(run_id, mesh_recipe), None


def _predict_next(elapsed_history: list[float], fallback_growth: float) -> float | None:
  if not elapsed_history:
    return None
  if len(elapsed_history) == 1:
    return elapsed_history[-1] * fallback_growth
  ratio = elapsed_history[-1] / max(elapsed_history[-2], 1.0e-12)
  return elapsed_history[-1] * max(ratio, fallback_growth)


def _write_table(df: pd.DataFrame, path: Path) -> None:
  df.to_csv(path.with_suffix(".csv"), index=False)
  df.to_pickle(path.with_suffix(".pkl"))


def _add_sequence_diagnostics(summary: pd.DataFrame, growth: int) -> pd.DataFrame:
  """Add post-facto sequence diagnostics for spatial-refinement studies."""

  if summary.empty:
    return summary
  df = summary.copy()
  df["h"] = df["spacing"]
  df["delta_pcm"] = df["delta_keff"] * 1.0e5
  df["keff_delta_from_previous"] = np.nan
  df["keff_delta_ratio"] = np.nan
  df["apparent_order_abs"] = np.nan
  df["aitken_keff"] = np.nan
  df["aitken_delta_from_current"] = np.nan

  keys = [
    "case",
    "index",
    "polar",
    "azimuth",
    "quadrature",
    "spacing_mode",
    "density_correction",
  ]
  completed = df["status"] == "completed"
  for _, group in df[completed].sort_values("mesh_n").groupby(keys, dropna=False):
    indices = list(group.index)
    values = [float(df.loc[index, "keff"]) for index in indices]
    for local, index in enumerate(indices):
      if local >= 1:
        delta = values[local] - values[local - 1]
        df.loc[index, "keff_delta_from_previous"] = delta
      if local >= 2:
        previous_delta = values[local - 1] - values[local - 2]
        current_delta = values[local] - values[local - 1]
        if abs(previous_delta) > 0.0:
          ratio = current_delta / previous_delta
          df.loc[index, "keff_delta_ratio"] = ratio
          if abs(ratio) > 0.0:
            df.loc[index, "apparent_order_abs"] = -math.log(abs(ratio), growth)
        denom = values[local] - 2.0 * values[local - 1] + values[local - 2]
        if abs(denom) > 0.0:
          aitken = values[local] - (values[local] - values[local - 1]) ** 2 / denom
          df.loc[index, "aitken_keff"] = aitken
          df.loc[index, "aitken_delta_from_current"] = aitken - values[local]
  return df


def run_sweep(args) -> dict:
  timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
  output_dir = Path(args.output_dir) if args.output_dir else (
    Path(__file__).resolve().parent / "results" / "moc_mesh_sweeps" / timestamp
  )
  output_dir.mkdir(parents=True, exist_ok=True)

  manifest = {
    "created_utc": timestamp,
    "purpose": "C5G7 centered-mesh MOC pin-cell spatial doubling sweep",
    "case": args.case,
    "index": args.index,
    "mesh_start": args.mesh_start,
    "mesh_growth": args.mesh_growth,
    "mesh_max": args.mesh_max,
    "polar_values": args.polar_values,
    "azimuth_values": args.azimuth_values,
    "quadrature": args.quadrature,
    "spacing_mode": args.spacing_mode,
    "spacing": args.spacing,
    "spacing_scale": args.spacing_scale,
    "density_correction": args.density_correction,
    "predicted_stop_seconds": args.predicted_stop_seconds,
    "threads": args.threads,
    "environment": {
      "python": sys.version,
      "executable": sys.executable,
      "platform": platform.platform(),
      "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
      "pythonpath": os.environ.get("PYTHONPATH"),
      "ld_library_path": os.environ.get("LD_LIBRARY_PATH"),
      "git_commit": _git_commit(),
    },
    "solver": {
      "equation": args.equation,
      "inner_iters": args.inner_iters,
      "outer_iters": args.outer_iters,
      "eigen_iters": args.eigen_iters,
      "inner_tolerance": args.inner_tolerance,
      "outer_tolerance": args.outer_tolerance,
      "eigen_tolerance": args.eigen_tolerance,
      "model": "flattened: one inner iteration and one outer iteration per k-update by default",
    },
  }
  (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

  summary_rows: list[dict] = []
  cell_flux_tables: list[pd.DataFrame] = []
  material_flux_tables: list[pd.DataFrame] = []
  material_area_rows: list[dict] = []

  for polar in args.polar_values:
    for azimuth in args.azimuth_values:
      elapsed_history: list[float] = []
      for mesh_n in _mesh_values(args.mesh_start, args.mesh_growth, args.mesh_max):
        predicted = _predict_next(
          elapsed_history,
          fallback_growth=args.mesh_growth * args.mesh_growth,
        )
        if predicted is not None and predicted > args.predicted_stop_seconds:
          stop = SweepRun(
            run_id=f"pin{args.index}_n{mesh_n}_p{polar}_a{azimuth}_skipped",
            case=args.case,
            index=args.index,
            mesh_n=mesh_n,
            cells=mesh_n * mesh_n,
            azimuth=azimuth,
            polar=polar,
            quadrature=args.quadrature,
            spacing=float("nan"),
            spacing_mode=args.spacing_mode,
            density_correction=bool(args.density_correction),
            equation=args.equation,
            elapsed_s=None,
            predicted_next_s=predicted,
            stop_reason="predicted_runtime_exceeds_limit",
            status="skipped",
          )
          summary_rows.append(asdict(stop))
          break

        print(
          f"[sweep] running mesh={mesh_n} polar={polar} azimuth={azimuth}",
          flush=True,
        )
        run, cell_df, material_df, area_rows, elapsed = _run_case(
          args,
          mesh_n,
          polar,
          azimuth,
          output_dir,
        )
        predicted_next = _predict_next(
          elapsed_history + ([elapsed] if elapsed is not None else []),
          fallback_growth=args.mesh_growth * args.mesh_growth,
        )
        row = asdict(run)
        row["predicted_next_s"] = predicted_next
        summary_rows.append(row)
        if not cell_df.empty:
          cell_flux_tables.append(cell_df)
        if not material_df.empty:
          material_flux_tables.append(material_df)
        material_area_rows.extend(area_rows)
        if elapsed is None:
          break
        elapsed_history.append(elapsed)

  summary = _add_sequence_diagnostics(pd.DataFrame(summary_rows), args.mesh_growth)
  cell_flux = pd.concat(cell_flux_tables, ignore_index=True) if cell_flux_tables else pd.DataFrame()
  material_flux = (
    pd.concat(material_flux_tables, ignore_index=True)
    if material_flux_tables else pd.DataFrame()
  )
  material_area = pd.DataFrame(material_area_rows)

  _write_table(summary, output_dir / "summary")
  _write_table(cell_flux, output_dir / "cell_flux")
  _write_table(material_flux, output_dir / "material_flux")
  _write_table(material_area, output_dir / "material_area_corrections")

  artifacts = {
    "output_dir": str(output_dir),
    "manifest": str(output_dir / "manifest.json"),
    "summary_csv": str(output_dir / "summary.csv"),
    "summary_pkl": str(output_dir / "summary.pkl"),
    "cell_flux_csv": str(output_dir / "cell_flux.csv"),
    "cell_flux_pkl": str(output_dir / "cell_flux.pkl"),
    "material_flux_csv": str(output_dir / "material_flux.csv"),
    "material_flux_pkl": str(output_dir / "material_flux.pkl"),
    "material_area_corrections_csv": str(output_dir / "material_area_corrections.csv"),
    "material_area_corrections_pkl": str(output_dir / "material_area_corrections.pkl"),
  }
  (output_dir / "artifacts.json").write_text(json.dumps(artifacts, indent=2, sort_keys=True) + "\n")
  return artifacts


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
  parser.add_argument("--case", choices=("pin",), default="pin")
  parser.add_argument("--index", type=int, default=0)
  parser.add_argument("--mesh-start", type=int, default=4)
  parser.add_argument("--mesh-growth", type=int, default=2)
  parser.add_argument("--mesh-max", type=int, default=4096)
  parser.add_argument("--polar-values", type=_parse_int_list, default=[2])
  parser.add_argument("--azimuth-values", type=_parse_int_list, default=[4])
  parser.add_argument("--quadrature", default="u-gl")
  parser.add_argument("--spacing-mode", choices=("cell", "uniform"), default="cell")
  parser.add_argument("--spacing", type=float, default=0.25)
  parser.add_argument("--spacing-scale", type=float, default=1.0)
  parser.add_argument("--density-correction", dest="density_correction",
                      action="store_true", default=True)
  parser.add_argument("--no-density-correction", dest="density_correction",
                      action="store_false")
  parser.add_argument("--predicted-stop-seconds", type=float, default=600.0)
  parser.add_argument("--threads", type=int, default=16)
  parser.add_argument("--equation", choices=("scmoc", "ddmoc"), default="scmoc")
  parser.add_argument("--inner-tolerance", type=float, default=1.0e-12)
  parser.add_argument("--outer-tolerance", type=float, default=1.0e-12)
  parser.add_argument("--eigen-tolerance", type=float, default=1.0e-5)
  parser.add_argument("--inner-iters", type=int, default=1)
  parser.add_argument("--outer-iters", type=int, default=1)
  parser.add_argument("--eigen-iters", type=int, default=2000)
  parser.add_argument("--monitor-level", type=int, default=0)
  parser.add_argument("--print-interval", type=int, default=25)
  parser.add_argument("--output-dir", default=None)
  return parser


def main() -> None:
  parser = add_arguments(
    argparse.ArgumentParser(description="Run a C5G7 centered-mesh MOC doubling sweep.")
  )
  args = parser.parse_args()
  os.environ["OMP_NUM_THREADS"] = str(args.threads)
  artifacts = run_sweep(args)
  print(json.dumps(artifacts, indent=2, sort_keys=True))


if __name__ == "__main__":
  main()
