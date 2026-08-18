"""Run true CSG-MOC overlay sweeps for the C5G7 UO2 pin cell."""

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
  from .csg_c5g7 import (
    PinDiscretization,
    build_pin,
    run_native_csg_moc,
  )
except ImportError:
  from csg_c5g7 import (
    PinDiscretization,
    build_pin,
    run_native_csg_moc,
  )


@dataclass(frozen=True)
class CSGMOCSweepRun:
  run_id: str
  case: str
  index: int
  mesh_n: int
  cells: int | None
  fsrs: int | None
  azimuth: int
  polar: int
  quadrature: str
  spacing: float
  overlay_variant: str
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


def _parse_float_list(spec: str) -> list[float]:
  values = []
  for part in str(spec).split(","):
    part = part.strip()
    if part:
      values.append(float(part))
  if not values:
    raise argparse.ArgumentTypeError("list must contain at least one value")
  return values


def _mesh_values(start: int, growth: int, max_mesh: int) -> Iterable[int]:
  current = int(start)
  while current <= int(max_mesh):
    yield current
    current *= int(growth)


def _git_commit() -> str | None:
  try:
    return subprocess.check_output(
      ["git", "rev-parse", "HEAD"],
      text=True,
      stderr=subprocess.DEVNULL,
    ).strip()
  except Exception:
    return None


def _predict_next(elapsed_history: list[float], fallback_growth: float) -> float | None:
  if not elapsed_history:
    return None
  if len(elapsed_history) == 1:
    return elapsed_history[-1] * fallback_growth
  ratio = elapsed_history[-1] / max(elapsed_history[-2], 1.0e-12)
  return elapsed_history[-1] * max(ratio, fallback_growth)


def _run_args(args, mesh_n: int, polar: int, azimuth: int, spacing: float):
  return argparse.Namespace(
    csg_overlay_mesh_variant=args.overlay_variant,
    csg_overlay_pin_mesh=mesh_n,
    moc_equation=args.equation,
    moc_polar=polar,
    moc_azimuth=azimuth,
    moc_quadrature=args.quadrature,
    moc_spacing=spacing,
    moc_inner_tolerance=args.inner_tolerance,
    moc_outer_tolerance=args.outer_tolerance,
    moc_eigen_tolerance=args.eigen_tolerance,
    moc_inner_iters=args.inner_iters,
    moc_outer_iters=args.outer_iters,
    moc_eigen_iters=args.eigen_iters,
    moc_monitor_level=args.monitor_level,
    moc_print_interval=args.print_interval,
  )


def _run_case(args, csg_case, mesh_n: int, polar: int, azimuth: int, spacing: float, output_dir: Path):
  spacing_label = f"{spacing:.8g}".replace(".", "p").replace("-", "m")
  run_id = f"csg_moc_pin{args.index}_n{mesh_n}_p{polar}_a{azimuth}_s{spacing_label}"
  run_dir = output_dir / "runs"
  run_dir.mkdir(parents=True, exist_ok=True)
  result_path = run_dir / f"{run_id}.json"
  run_args = _run_args(args, mesh_n, polar, azimuth, spacing)

  start = time.time()
  try:
    transport = run_native_csg_moc(csg_case, run_args)
    elapsed = time.time() - start
    transport.update({
      "run_id": run_id,
      "status": "completed",
      "elapsed_s": elapsed,
    })
    result_path.write_text(json.dumps(transport, indent=2, sort_keys=True) + "\n")
    overlay_product = transport.get("overlay_product") or {}
    run = CSGMOCSweepRun(
      run_id=run_id,
      case=args.case,
      index=args.index,
      mesh_n=mesh_n,
      cells=int(transport.get("cells")) if transport.get("cells") is not None else None,
      fsrs=int(overlay_product.get("fsrs")) if overlay_product.get("fsrs") is not None else None,
      azimuth=azimuth,
      polar=polar,
      quadrature=args.quadrature,
      spacing=spacing,
      overlay_variant=args.overlay_variant,
      equation=args.equation,
      elapsed_s=elapsed,
      predicted_next_s=None,
      stop_reason=None,
      status="completed",
      keff=float(transport["keff"]),
      finite=bool(transport["finite"]),
      sweeps=int(transport.get("sweeps", 0)),
      reference_keff=(
        float(transport["reference_keff"])
        if transport.get("reference_keff") is not None else None
      ),
      delta_keff=(
        float(transport["delta_keff"])
        if transport.get("delta_keff") is not None else None
      ),
      result_path=str(result_path),
    )
    return run, elapsed
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
      "spacing": spacing,
      "overlay_variant": args.overlay_variant,
      "equation": args.equation,
      "elapsed_s": elapsed,
    }
    result_path.write_text(json.dumps(failure, indent=2, sort_keys=True) + "\n")
    run = CSGMOCSweepRun(
      run_id=run_id,
      case=args.case,
      index=args.index,
      mesh_n=mesh_n,
      cells=None,
      fsrs=None,
      azimuth=azimuth,
      polar=polar,
      quadrature=args.quadrature,
      spacing=spacing,
      overlay_variant=args.overlay_variant,
      equation=args.equation,
      elapsed_s=elapsed,
      predicted_next_s=None,
      stop_reason="failed",
      status="failed",
      error=repr(exc),
      result_path=str(result_path),
    )
    return run, None


def _write_table(df: pd.DataFrame, path: Path) -> None:
  df.to_csv(path.with_suffix(".csv"), index=False)
  df.to_pickle(path.with_suffix(".pkl"))


def _add_sequence_diagnostics(summary: pd.DataFrame, growth: int) -> pd.DataFrame:
  if summary.empty:
    return summary
  df = summary.copy()
  df["h"] = 1.0 / df["mesh_n"].astype(float)
  df["delta_pcm"] = df["delta_keff"] * 1.0e5
  df["keff_delta_from_previous"] = np.nan
  df["keff_delta_ratio"] = np.nan
  df["apparent_order_abs"] = np.nan
  df["aitken_keff"] = np.nan
  df["aitken_delta_from_current"] = np.nan
  keys = ["case", "index", "polar", "azimuth", "quadrature", "spacing", "overlay_variant"]
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
    Path(__file__).resolve().parent / "results" / "csg_moc_overlay_sweeps" / timestamp
  )
  output_dir.mkdir(parents=True, exist_ok=True)

  discretization = PinDiscretization(
    fuel_radial_divisions=args.fuel_radial_divisions,
    moderator_radial_divisions=args.moderator_radial_divisions,
    azimuthal_slices=args.azimuthal_slices,
    azimuthal_offset_degrees=args.azimuthal_offset_degrees,
  )
  csg_case = build_pin(args.index, discretization)

  manifest = {
    "created_utc": timestamp,
    "purpose": "C5G7 UO2 pin true CSG-MOC centered-overlay sweep",
    "case": args.case,
    "index": args.index,
    "mesh_start": args.mesh_start,
    "mesh_growth": args.mesh_growth,
    "mesh_max": args.mesh_max,
    "polar_values": args.polar_values,
    "azimuth_values": args.azimuth_values,
    "spacing_values": args.spacing_values,
    "quadrature": args.quadrature,
    "overlay_variant": args.overlay_variant,
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
    "csg_discretization": asdict(discretization),
    "solver": {
      "equation": args.equation,
      "inner_iters": args.inner_iters,
      "outer_iters": args.outer_iters,
      "eigen_iters": args.eigen_iters,
      "inner_tolerance": args.inner_tolerance,
      "outer_tolerance": args.outer_tolerance,
      "eigen_tolerance": args.eigen_tolerance,
      "model": "flattened power iteration when inner_iters=outer_iters=1",
    },
  }
  (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

  summary_rows: list[dict] = []
  for polar in args.polar_values:
    for azimuth in args.azimuth_values:
      for spacing in args.spacing_values:
        elapsed_history: list[float] = []
        for mesh_n in _mesh_values(args.mesh_start, args.mesh_growth, args.mesh_max):
          predicted = _predict_next(
            elapsed_history,
            fallback_growth=args.mesh_growth * args.mesh_growth,
          )
          if predicted is not None and predicted > args.predicted_stop_seconds:
            stop = CSGMOCSweepRun(
              run_id=f"csg_moc_pin{args.index}_n{mesh_n}_p{polar}_a{azimuth}_skipped",
              case=args.case,
              index=args.index,
              mesh_n=mesh_n,
              cells=None,
              fsrs=None,
              azimuth=azimuth,
              polar=polar,
              quadrature=args.quadrature,
              spacing=spacing,
              overlay_variant=args.overlay_variant,
              equation=args.equation,
              elapsed_s=None,
              predicted_next_s=predicted,
              stop_reason="predicted_runtime_exceeds_limit",
              status="skipped",
            )
            summary_rows.append(asdict(stop))
            break
          print(
            f"[sweep] running csg-moc mesh={mesh_n} polar={polar} azimuth={azimuth} spacing={spacing}",
            flush=True,
          )
          run, elapsed = _run_case(args, csg_case, mesh_n, polar, azimuth, spacing, output_dir)
          predicted_next = _predict_next(
            elapsed_history + ([elapsed] if elapsed is not None else []),
            fallback_growth=args.mesh_growth * args.mesh_growth,
          )
          row = asdict(run)
          row["predicted_next_s"] = predicted_next
          summary_rows.append(row)
          if elapsed is None:
            break
          elapsed_history.append(elapsed)

  summary = _add_sequence_diagnostics(pd.DataFrame(summary_rows), args.mesh_growth)
  _write_table(summary, output_dir / "summary")
  artifacts = {
    "output_dir": str(output_dir),
    "manifest": str(output_dir / "manifest.json"),
    "summary_csv": str(output_dir / "summary.csv"),
    "summary_pkl": str(output_dir / "summary.pkl"),
  }
  (output_dir / "artifacts.json").write_text(json.dumps(artifacts, indent=2, sort_keys=True) + "\n")
  return artifacts


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
  parser.add_argument("--case", choices=("pin",), default="pin")
  parser.add_argument("--index", type=int, default=0)
  parser.add_argument("--mesh-start", type=int, default=4)
  parser.add_argument("--mesh-growth", type=int, default=2)
  parser.add_argument("--mesh-max", type=int, default=512)
  parser.add_argument("--polar-values", type=_parse_int_list, default=[2])
  parser.add_argument("--azimuth-values", type=_parse_int_list, default=[4])
  parser.add_argument("--spacing-values", type=_parse_float_list, default=[0.02])
  parser.add_argument("--quadrature", default="u-gl")
  parser.add_argument("--overlay-variant", choices=("legacy", "optimized", "centered"), default="centered")
  parser.add_argument("--predicted-stop-seconds", type=float, default=600.0)
  parser.add_argument("--threads", type=int, default=16)
  parser.add_argument("--equation", choices=("scmoc", "ddmoc"), default="ddmoc")
  parser.add_argument("--inner-tolerance", type=float, default=1.0e-12)
  parser.add_argument("--outer-tolerance", type=float, default=1.0e-12)
  parser.add_argument("--eigen-tolerance", type=float, default=1.0e-8)
  parser.add_argument("--inner-iters", type=int, default=1)
  parser.add_argument("--outer-iters", type=int, default=1)
  parser.add_argument("--eigen-iters", type=int, default=8000)
  parser.add_argument("--monitor-level", type=int, default=0)
  parser.add_argument("--print-interval", type=int, default=50)
  parser.add_argument("--fuel-radial-divisions", type=int, default=1)
  parser.add_argument("--moderator-radial-divisions", type=int, default=1)
  parser.add_argument("--azimuthal-slices", type=int, default=8)
  parser.add_argument("--azimuthal-offset-degrees", type=float, default=0.0)
  parser.add_argument("--output-dir", default=None)
  return parser


def main() -> None:
  parser = add_arguments(
    argparse.ArgumentParser(description="Run a C5G7 UO2 pin true CSG-MOC overlay sweep.")
  )
  args = parser.parse_args()
  os.environ["OMP_NUM_THREADS"] = str(args.threads)
  artifacts = run_sweep(args)
  print(json.dumps(artifacts, indent=2, sort_keys=True))


if __name__ == "__main__":
  main()
