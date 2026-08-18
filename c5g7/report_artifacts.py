"""Generate report figures and pin-error rows for the C5G7 CSG report."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

try:
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt
  from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
  from matplotlib.patches import Patch
except ImportError as exc:
  raise RuntimeError("matplotlib is required for report artifact generation") from exc

try:
  from .csg_c5g7 import build_assembly, build_core, build_pin
except ImportError:
  from csg_c5g7 import build_assembly, build_core, build_pin


MATERIAL_COLORS = {
  0: "#2f6db5",
  1: "#c77c2e",
  2: "#b85827",
  3: "#8f2d20",
  4: "#7aa6c2",
  5: "#5b4a9b",
  6: "#d7dde5",
}

MATERIAL_LABELS = {
  0: "UO2",
  1: "MOX 4.3%",
  2: "MOX 7.0%",
  3: "MOX 8.7%",
  4: "Guide tube",
  5: "Fission chamber",
  6: "Moderator",
}


def _material_cmap():
  ids = sorted(MATERIAL_COLORS)
  colors = [MATERIAL_COLORS[i] for i in ids]
  cmap = ListedColormap(colors)
  boundaries = [i - 0.5 for i in ids] + [ids[-1] + 0.5]
  norm = BoundaryNorm(boundaries, cmap.N)
  return cmap, norm


def _plot_material_slice(ax, csg_case, nx: int, title: str) -> None:
  product = csg_case.model.sample_slice(nx, nx)
  cmap, norm = _material_cmap()
  ax.imshow(
    product.material_id,
    origin="lower",
    interpolation="nearest",
    cmap=cmap,
    norm=norm,
    extent=[0.0, csg_case.model.width_x, 0.0, csg_case.model.width_y],
  )
  ax.set_title(title)
  ax.set_xlabel("x [cm]")
  ax.set_ylabel("y [cm]")
  ax.set_aspect("equal")


def _hashed_region_field(region_id: np.ndarray) -> np.ndarray:
  ids = region_id.astype(np.uint64, copy=False)
  return ((ids * np.uint64(2654435761)) % np.uint64(2**32)).astype(float) / float(2**32)


def _region_boundary_mask(region_id: np.ndarray) -> np.ndarray:
  boundary = np.zeros(region_id.shape, dtype=bool)
  boundary[1:, :] |= region_id[1:, :] != region_id[:-1, :]
  boundary[:-1, :] |= region_id[:-1, :] != region_id[1:, :]
  boundary[:, 1:] |= region_id[:, 1:] != region_id[:, :-1]
  boundary[:, :-1] |= region_id[:, :-1] != region_id[:, 1:]
  return boundary


def _plot_region_slice(
  ax,
  csg_case,
  nx: int,
  title: str,
  *,
  xlim: tuple[float, float] | None = None,
  ylim: tuple[float, float] | None = None,
) -> None:
  product = csg_case.model.sample_slice(nx, nx)
  region_id = np.asarray(product.region_id)
  extent = [0.0, csg_case.model.width_x, 0.0, csg_case.model.width_y]
  ax.imshow(
    _hashed_region_field(region_id),
    origin="lower",
    interpolation="nearest",
    cmap="turbo",
    extent=extent,
  )
  boundary = np.ma.masked_where(~_region_boundary_mask(region_id), np.ones(region_id.shape))
  ax.imshow(
    boundary,
    origin="lower",
    interpolation="nearest",
    cmap=ListedColormap(["#111111"]),
    alpha=0.82,
    extent=extent,
  )
  ax.set_title(title)
  ax.set_xlabel("x [cm]")
  ax.set_ylabel("y [cm]")
  ax.set_aspect("equal")
  if xlim is not None:
    ax.set_xlim(*xlim)
  if ylim is not None:
    ax.set_ylim(*ylim)


def _draw_pin_grid_cm(ax, xlim: tuple[float, float], ylim: tuple[float, float]) -> None:
  x0 = int(np.floor(xlim[0] / 1.26))
  x1 = int(np.ceil(xlim[1] / 1.26))
  y0 = int(np.floor(ylim[0] / 1.26))
  y1 = int(np.ceil(ylim[1] / 1.26))
  for i in range(x0, x1 + 1):
    ax.axvline(i * 1.26, color="white", linewidth=0.4, alpha=0.8)
  for j in range(y0, y1 + 1):
    ax.axhline(j * 1.26, color="white", linewidth=0.4, alpha=0.8)


def write_geometry_figure(output_dir: Path) -> Path:
  """Write Detran-sampled CSG material maps for pin, assembly, and core."""

  output_dir.mkdir(parents=True, exist_ok=True)
  path = output_dir / "csg_geometry_material_maps.pdf"
  cases = [
    (build_pin(0), 180, "UO2 pin"),
    (build_assembly(0), 680, "UO2 assembly"),
    (build_core(), 612, "Quarter core"),
  ]

  with matplotlib.rc_context({"font.size": 8}):
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65))
    for ax, (case, samples, title) in zip(axes, cases):
      _plot_material_slice(ax, case, samples, title)
    legend_handles = [
      Patch(facecolor=MATERIAL_COLORS[i], edgecolor="none", label=MATERIAL_LABELS[i])
      for i in sorted(MATERIAL_COLORS)
    ]
    fig.legend(
      handles=legend_handles,
      loc="upper center",
      ncol=4,
      frameon=False,
      bbox_to_anchor=(0.5, 1.03),
    )
    fig.subplots_adjust(top=0.72, wspace=0.28)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
  return path


def write_flat_source_region_figure(output_dir: Path) -> Path:
  """Write Detran-sampled flat-source region maps for representative CSG cases."""

  output_dir.mkdir(parents=True, exist_ok=True)
  path = output_dir / "csg_flat_source_regions.pdf"
  pin_pitch = 1.26
  central_patch = (6 * pin_pitch, 11 * pin_pitch)
  cases = [
    (build_pin(0), 360, "UO2 pin", None, None),
    (build_assembly(0), 1020, "UO2 assembly FSRs", central_patch, central_patch),
    (build_assembly(1), 1020, "MOX assembly FSRs", central_patch, central_patch),
  ]

  with matplotlib.rc_context({"font.size": 8}):
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), constrained_layout=True)
    for ax, (case, samples, title, xlim, ylim) in zip(axes, cases):
      _plot_region_slice(ax, case, samples, title, xlim=xlim, ylim=ylim)
      if xlim is not None and ylim is not None:
        _draw_pin_grid_cm(ax, xlim, ylim)
      ax.text(
        0.02,
        0.02,
        f"{case.model.number_regions:,} regions",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7,
        color="white",
        bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.58, "pad": 2},
      )
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
  return path


def _masked_reference(array: np.ndarray, reference: np.ndarray) -> np.ma.MaskedArray:
  return np.ma.masked_where(reference <= 0.0, array)


def _draw_pin_grid(ax, pins_per_side: int = 51) -> None:
  for value in range(0, pins_per_side + 1, 17):
    ax.axhline(value - 0.5, color="black", linewidth=0.55, alpha=0.8)
    ax.axvline(value - 0.5, color="black", linewidth=0.55, alpha=0.8)
  ax.set_xticks([8, 25, 42])
  ax.set_yticks([8, 25, 42])
  ax.set_xticklabels(["A0", "A1", "Refl."])
  ax.set_yticklabels(["A0", "A1", "Refl."])
  ax.set_xlim(-0.5, pins_per_side - 0.5)
  ax.set_ylim(-0.5, pins_per_side - 0.5)


def write_pin_power_figure(core_result_path: Path, output_dir: Path) -> tuple[Path, dict]:
  """Write computed/reference/error core pin-power figure from a CSG result."""

  output_dir.mkdir(parents=True, exist_ok=True)
  result = json.loads(core_result_path.read_text())
  pin_power = result["transport"]["pin_power"]
  calculated = np.asarray(pin_power["calculated_normalized"], dtype=float)
  reference = np.asarray(pin_power["reference"], dtype=float)
  relative_error = np.asarray(pin_power["relative_error_percent"], dtype=float)

  path = output_dir / "csg_core_pin_power_errors.pdf"
  max_power = float(max(np.nanmax(calculated), np.nanmax(reference)))
  max_abs_error = float(np.nanmax(np.abs(relative_error)))
  max_abs_error = max(max_abs_error, 1.0)

  with matplotlib.rc_context({"font.size": 8}):
    fig, axes = plt.subplots(1, 3, figsize=(7.25, 2.45), constrained_layout=True)
    panels = [
      (
        _masked_reference(reference, reference),
        "Reference power",
        "magma",
        None,
        0.0,
        max_power,
      ),
      (
        _masked_reference(calculated, reference),
        "CSG-MOC power",
        "magma",
        None,
        0.0,
        max_power,
      ),
      (
        _masked_reference(relative_error, reference),
        "Relative error [%]",
        "coolwarm",
        TwoSlopeNorm(vmin=-max_abs_error, vcenter=0.0, vmax=max_abs_error),
        None,
        None,
      ),
    ]
    for ax, (data, title, cmap, norm, vmin, vmax) in zip(axes, panels):
      image = ax.imshow(
        data,
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
      )
      ax.set_title(title)
      ax.set_aspect("equal")
      _draw_pin_grid(ax, reference.shape[0])
      fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
  return path, pin_power


def _format_float(value, digits=4) -> str:
  return f"{float(value):.{digits}g}"


def _latest_summary(root: Path) -> Path | None:
  paths = sorted(root.glob("*/summary.csv"))
  return paths[-1] if paths else None


def _read_completed_sweep_rows(path: Path, method: str) -> list[dict]:
  rows = []
  if path is None or not path.exists():
    return rows
  with path.open(newline="") as stream:
    for row in csv.DictReader(stream):
      if row.get("status") != "completed":
        continue
      if row.get("case") != "pin" or row.get("index") != "0":
        continue
      row = dict(row)
      row["method"] = method
      rows.append(row)
  return rows


def _format_table_float(value: str, digits: int = 4) -> str:
  if value is None or value == "":
    return "--"
  return f"{float(value):.{digits}g}"


def _format_keff(value: str) -> str:
  if value is None or value == "":
    return "--"
  return f"{float(value):.6f}"


def _format_pcm(value: str) -> str:
  if value is None or value == "":
    return "--"
  return f"{float(value):.1f}"


def _format_int(value: str) -> str:
  if value is None or value == "":
    return "--"
  return f"{int(float(value)):,}"


def _sweep_sort_key(row: dict) -> tuple[int, int, int]:
  method_order = {"MOC": 0, "SN": 1}.get(row["method"], 9)
  equation_order = {"scmoc": 0, "dd": 1, "sc": 2}.get(row.get("equation", ""), 9)
  return method_order, equation_order, int(row["mesh_n"])


def _sweep_equation_label(row: dict) -> str:
  equation = row.get("equation", "")
  if equation == "scmoc":
    return "SC-MOC"
  return equation.upper()


def write_uo2_pin_mesh_sweep_table(
  output_dir: Path,
  moc_summary: Path | None,
  sn_summary: Path | None,
) -> Path:
  """Write a LaTeX table for the NxN UO2 pin-cell MOC/SN sweeps."""

  output_dir.mkdir(parents=True, exist_ok=True)
  table_path = output_dir / "mesh_sweep_uo2_pin_table.tex"
  rows = _read_completed_sweep_rows(moc_summary, "MOC")
  rows.extend(_read_completed_sweep_rows(sn_summary, "SN"))
  lines = []
  if moc_summary is not None:
    lines.append(f"% MOC summary: {moc_summary}")
  if sn_summary is not None:
    lines.append(f"% SN summary: {sn_summary}")
  lines.extend([
    r"\begin{center}",
    r"\tiny",
    r"\setlength{\tabcolsep}{2.5pt}",
    r"\begin{tabular}{llrrrrrrr}",
    r"\toprule",
    r"Method & Eq. & \(N\) & Cells & \(h\) [cm] & Detran \(k\) & pcm & Time [s] & Sweeps \\",
    r"\midrule",
  ])
  for row in sorted(rows, key=_sweep_sort_key):
    h_value = row.get("spacing") or row.get("cell_width")
    lines.append(" & ".join([
      row["method"],
      _sweep_equation_label(row),
      _format_int(row["mesh_n"]),
      _format_int(row["cells"]),
      _format_table_float(h_value, 5),
      _format_keff(row.get("keff")),
      _format_pcm(row.get("delta_pcm")),
      _format_table_float(row.get("elapsed_s"), 4),
      _format_int(row.get("sweeps")),
    ]) + r"\\")
  lines.extend([
    r"\bottomrule",
    r"\end{tabular}",
    r"\end{center}",
  ])
  table_path.write_text("\n".join(lines) + "\n")
  return table_path


def _read_json(path: Path) -> dict | None:
  if path is None or not path.exists():
    return None
  return json.loads(path.read_text())


def _moc_accuracy_row(path: Path) -> list[str] | None:
  result = _read_json(path)
  if result is None:
    return None
  recipe = result.get("pin_mesh_recipe", {})
  mesh = "{0}x{0} {1}".format(
    recipe.get("pin_mesh", 7),
    "legacy",
  )
  delta_pcm = 1.0e5 * float(result["delta_keff"])
  return [
    "MOC",
    mesh,
    result["equation"].upper(),
    "asqr-asdr",
    str(result["polar"]),
    str(result["azimuth"]),
    _format_table_float(result.get("spacing"), 4),
    _format_keff(result.get("keff")),
    _format_pcm(delta_pcm),
    _format_table_float(result.get("elapsed"), 4),
    _format_int(result.get("sweeps")),
  ]


def _sn_accuracy_row(path: Path) -> list[str] | None:
  result = _read_json(path)
  if result is None:
    return None
  recipe = result.get("pin_mesh_recipe", {})
  discretization = result["discretization"]
  values = result["results"]
  mesh = "{0}x{0} optimized".format(recipe.get("number", 7))
  return [
    "SN",
    mesh,
    discretization["equation"].upper(),
    discretization["quadrature"],
    str(discretization["polar_per_octant"]),
    str(discretization["azimuth_per_octant"]),
    "--",
    _format_keff(values.get("keff")),
    _format_pcm(values.get("delta_pcm")),
    _format_table_float(values.get("elapsed_s"), 4),
    _format_int(values.get("sweeps")),
  ]


def write_uo2_pin_accuracy_table(
  output_dir: Path,
  moc_result: Path | None,
  sn_result: Path | None,
) -> Path:
  """Write the few-pcm UO2 pin-cell MOC/SN accuracy table."""

  output_dir.mkdir(parents=True, exist_ok=True)
  table_path = output_dir / "uo2_pin_accuracy_table.tex"
  rows = []
  moc_row = _moc_accuracy_row(moc_result)
  if moc_row is not None:
    rows.append(moc_row)
  sn_row = _sn_accuracy_row(sn_result)
  if sn_row is not None:
    rows.append(sn_row)
  lines = []
  if moc_result is not None:
    lines.append(f"% MOC result: {moc_result}")
  if sn_result is not None:
    lines.append(f"% SN result: {sn_result}")
  lines.extend([
    r"\begin{center}",
    r"\scriptsize",
    r"\setlength{\tabcolsep}{3pt}",
    r"\begin{tabular}{lllrrrrrrrr}",
    r"\toprule",
    r"Method & Mesh & Eq. & Quad. & Pol. & Az. & Spacing [cm] & Detran \(k\) & pcm & Time [s] & Sweeps \\",
    r"\midrule",
  ])
  for row in rows:
    lines.append(" & ".join(row) + r"\\")
  lines.extend([
    r"\bottomrule",
    r"\end{tabular}",
    r"\end{center}",
  ])
  table_path.write_text("\n".join(lines) + "\n")
  return table_path


def write_pin_error_row(pin_power: dict, output_dir: Path) -> Path:
  output_dir.mkdir(parents=True, exist_ok=True)
  row_path = output_dir / "csg_core_0_pin_error_row.tex"
  max_pin = pin_power["max_abs_pin"]
  row = " & ".join([
    "CSG-MOC",
    str(pin_power["active_pin_count"]),
    _format_float(pin_power["rms_percent"]),
    _format_float(pin_power["mean_abs_percent"]),
    _format_float(pin_power["max_abs_percent"]),
    f"({max_pin['i']}, {max_pin['j']})",
    _format_float(max_pin["reference"]),
    _format_float(max_pin["calculated"]),
  ])
  row_path.write_text(row + "\n")
  return row_path


def main(argv: list[str] | None = None) -> None:
  parser = argparse.ArgumentParser(
    description="Generate C5G7 report figures from Detran CSG artifacts."
  )
  parser.add_argument(
    "--csg-result-dir",
    type=Path,
    default=Path(__file__).resolve().parent / "csg_results",
  )
  parser.add_argument(
    "--figure-dir",
    type=Path,
    default=Path(__file__).resolve().parent / "results" / "figures",
  )
  parser.add_argument(
    "--moc-sweep-summary",
    type=Path,
    default=None,
  )
  parser.add_argument(
    "--sn-sweep-summary",
    type=Path,
    default=None,
  )
  parser.add_argument(
    "--moc-accuracy-result",
    type=Path,
    default=None,
  )
  parser.add_argument(
    "--sn-accuracy-result",
    type=Path,
    default=None,
  )
  args = parser.parse_args(argv)

  result_dir = Path(__file__).resolve().parent / "results"
  moc_summary = args.moc_sweep_summary or _latest_summary(result_dir / "moc_mesh_sweeps")
  sn_summary = args.sn_sweep_summary or _latest_summary(result_dir / "sn_mesh_sweeps")
  accuracy_dir = result_dir / "accuracy_probes"
  moc_accuracy = args.moc_accuracy_result or (
    accuracy_dir / "moc_legacy7_ddmoc_asqr_asdr_p4a8_s0p01.json"
  )
  sn_accuracy = args.sn_accuracy_result or (
    accuracy_dir / "sn_optimized_dd_asqr_asdr_p3a6.json"
  )

  geometry = write_geometry_figure(args.figure_dir)
  flat_source_regions = write_flat_source_region_figure(args.figure_dir)
  pin_figure, pin_power = write_pin_power_figure(
    args.csg_result_dir / "csg_core_0.json",
    args.figure_dir,
  )
  row = write_pin_error_row(pin_power, args.csg_result_dir)
  sweep_table = write_uo2_pin_mesh_sweep_table(result_dir, moc_summary, sn_summary)
  accuracy_table = write_uo2_pin_accuracy_table(result_dir, moc_accuracy, sn_accuracy)
  print(json.dumps({
    "geometry_figure": str(geometry),
    "flat_source_region_figure": str(flat_source_regions),
    "pin_power_figure": str(pin_figure),
    "pin_error_row": str(row),
    "uo2_pin_accuracy_table": str(accuracy_table),
    "uo2_pin_mesh_sweep_table": str(sweep_table),
    "pin_power_summary": {
      key: pin_power[key]
      for key in (
        "active_pin_count",
        "rms_percent",
        "mean_abs_percent",
        "max_abs_percent",
        "max_abs_pin",
      )
    },
  }, indent=2, sort_keys=True))


if __name__ == "__main__":
  main()
