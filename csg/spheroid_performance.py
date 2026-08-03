"""Unit-box spheroid CSG and tracking performance playground.

The cases here are intentionally small and analytic.  They are examples first,
but each step mirrors something a unit test can check: scale an axis-aligned
spheroid into a guarded unit box, slice it through a known plane, track it with
explicit directions, and compare track-derived material volumes to the exact
spheroid volume.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import time
from typing import Sequence

import numpy as np

import detran as dt


OUTER_MATERIAL = 0
SPHEROID_MATERIAL = 1
SPHEROID_REGION = 0
OUTER_REGION = 1

MATERIAL_NAMES = {
    OUTER_MATERIAL: "box_minus_spheroid",
    SPHEROID_MATERIAL: "spheroid",
}

DIRECTION_SETS = {
    "x": ((1.0, 0.0, 0.0),),
    "y": ((0.0, 1.0, 0.0),),
    "z": ((0.0, 0.0, 1.0),),
    "xyz": ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    "diagonal": ((1.0, 1.0, 1.0),),
}


def _triple(values: Sequence[float], name: str) -> tuple[float, float, float]:
    if isinstance(values, (str, bytes)):
        values = values.split(",")
    out = tuple(float(value) for value in values)
    if len(out) != 3:
        raise ValueError(f"{name} must have exactly three values")
    return out


def fit_spheroid_axes(
    center: Sequence[float],
    axis_ratios: Sequence[float],
    guard_min: float = 0.01,
    guard_max: float = 0.99,
) -> tuple[tuple[float, float, float], float]:
    """Scale axis ratios so the spheroid fits inside the guarded unit box."""

    center = _triple(center, "center")
    axis_ratios = _triple(axis_ratios, "axis_ratios")
    if guard_min < 0.0 or guard_max > 1.0 or guard_min >= guard_max:
        raise ValueError("guard bounds must satisfy 0 <= guard_min < guard_max <= 1")
    for value in center:
        if value <= guard_min or value >= guard_max:
            raise ValueError("spheroid center must be inside the guard bounds")
    for value in axis_ratios:
        if value <= 0.0:
            raise ValueError("spheroid axis ratios must be positive")

    margins = tuple(min(value - guard_min, guard_max - value) for value in center)
    scale = min(margin / axis for margin, axis in zip(margins, axis_ratios))
    return tuple(scale * axis for axis in axis_ratios), scale


@dataclass(frozen=True)
class SpheroidCase:
    """A hand-checkable spheroid embedded in the unit box."""

    name: str
    center: tuple[float, float, float] = (0.5, 0.5, 0.5)
    axis_ratios: tuple[float, float, float] = (1.0, 1.0, 1.0)
    guard_min: float = 0.01
    guard_max: float = 0.99
    spheroid_material: int = SPHEROID_MATERIAL
    outer_material: int = OUTER_MATERIAL

    @property
    def semi_axes(self) -> tuple[float, float, float]:
        axes, _ = fit_spheroid_axes(
            self.center, self.axis_ratios, self.guard_min, self.guard_max
        )
        return axes

    @property
    def scale(self) -> float:
        _, scale = fit_spheroid_axes(
            self.center, self.axis_ratios, self.guard_min, self.guard_max
        )
        return scale

    @property
    def extents(self) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
        return tuple(
            (center - axis, center + axis)
            for center, axis in zip(self.center, self.semi_axes)
        )

    @property
    def exact_volume(self) -> float:
        rx, ry, rz = self.semi_axes
        return 4.0 * math.pi * rx * ry * rz / 3.0

    def xy_slice_axes(self, z: float) -> tuple[float, float]:
        """Return the exact ellipse semi-axes in the x-y slice at ``z``."""

        rx, ry, rz = self.semi_axes
        dz = (float(z) - self.center[2]) / rz
        if abs(dz) > 1.0:
            return 0.0, 0.0
        factor = math.sqrt(max(0.0, 1.0 - dz * dz))
        return rx * factor, ry * factor

    def centerline_lengths(self, axis: int) -> tuple[float, float, float]:
        """Exact box/spheroid/box lengths through the spheroid center."""

        center = self.center[axis]
        radius = self.semi_axes[axis]
        return center - radius, 2.0 * radius, 1.0 - center - radius

    def model(self):
        x, y, z = self.center
        rx, ry, rz = self.semi_axes
        return (
            dt.CSGModel3D(1.0, 1.0, 1.0)
            .add_spheroid(self.spheroid_material, x, y, z, rx, ry, rz)
            .add_box_minus_spheroid(
                self.outer_material,
                0.0,
                1.0,
                0.0,
                1.0,
                0.0,
                1.0,
                x,
                y,
                z,
                rx,
                ry,
                rz,
            )
        )


CASES = {
    "centered_sphere": SpheroidCase("centered_sphere"),
    "oblong_x": SpheroidCase("oblong_x", axis_ratios=(1.0, 0.5, 0.25)),
    "offcenter_x": SpheroidCase(
        "offcenter_x", center=(0.35, 0.5, 0.5), axis_ratios=(1.0, 0.5, 0.25)
    ),
    "offcenter_y": SpheroidCase(
        "offcenter_y", center=(0.5, 0.65, 0.5), axis_ratios=(0.5, 1.0, 0.25)
    ),
    "offcenter_z": SpheroidCase(
        "offcenter_z", center=(0.5, 0.5, 0.35), axis_ratios=(0.5, 0.25, 1.0)
    ),
}


def estimate_region_volumes(product) -> dict[int, float]:
    """Estimate region volumes by length times track spatial weights."""

    direction_weights = product.directions.weight
    track_direction = product.tracks.direction_index
    track_weight = product.tracks.spatial_weight
    segment_track = product.segments.track_index
    segment_weight = (
        track_weight[segment_track] * direction_weights[track_direction[segment_track]]
    )
    total_direction_weight = float(np.sum(direction_weights))
    volumes = {}
    for region in np.unique(product.segments.region_id):
        mask = product.segments.region_id == region
        volumes[int(region)] = float(
            np.sum(product.segments.length[mask] * segment_weight[mask])
            / total_direction_weight
        )
    return volumes


def benchmark_case(
    case: SpheroidCase,
    *,
    directions: Sequence[Sequence[float]] = DIRECTION_SETS["xyz"],
    spacing: float = 0.05,
    minimum_tracks: int = 1,
    slice_cells: int = 80,
    repeats: int = 1,
    output_dir: str | Path | None = None,
) -> dict:
    """Build, slice, and track one case, returning timings and summary data."""

    timings = []
    model = None
    slice_product = None
    track_product = None
    for _ in range(int(repeats)):
        t0 = time.perf_counter()
        model = case.model()
        t1 = time.perf_counter()
        slice_product = model.sample_slice_xy(case.center[2], slice_cells, slice_cells)
        t2 = time.perf_counter()
        track_product = model.track(
            vectors=directions, spacing=spacing, minimum_tracks=minimum_tracks
        )
        t3 = time.perf_counter()
        timings.append(
            {
                "build_s": t1 - t0,
                "slice_s": t2 - t1,
                "track_s": t3 - t2,
                "total_s": t3 - t0,
            }
        )

    assert model is not None
    assert slice_product is not None
    assert track_product is not None

    volumes = estimate_region_volumes(track_product)
    spheroid_volume = volumes.get(SPHEROID_REGION, 0.0)
    total_tracked_volume = sum(volumes.values())
    analytic_slice_xy = model.quadratic_slices("xy", case.center[2])
    analytic_projection_xy = model.quadratic_projections("xy")
    result = {
        "name": case.name,
        "center": list(case.center),
        "axis_ratios": list(case.axis_ratios),
        "semi_axes": list(case.semi_axes),
        "scale": case.scale,
        "extents": [list(pair) for pair in case.extents],
        "slice_z": case.center[2],
        "slice_xy_axes": list(case.xy_slice_axes(case.center[2])),
        "analytic_slice_xy_bounds": (
            list(analytic_slice_xy[0].bounds) if analytic_slice_xy else None
        ),
        "analytic_projection_xy_bounds": (
            list(analytic_projection_xy[0].bounds) if analytic_projection_xy else None
        ),
        "slice_cells": int(slice_cells),
        "spacing": float(spacing),
        "minimum_tracks": int(minimum_tracks),
        "directions": [list(vector) for vector in directions],
        "tracks": int(len(track_product.tracks)),
        "segments": int(len(track_product.segments)),
        "exact_spheroid_volume": case.exact_volume,
        "exact_total_volume": 1.0,
        "tracked_region_volumes": volumes,
        "spheroid_volume_error": spheroid_volume - case.exact_volume,
        "total_volume_error": total_tracked_volume - 1.0,
        "timings": {
            key: float(np.mean([timing[key] for timing in timings]))
            for key in timings[0]
        },
    }

    if output_dir is not None:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        (output_path / f"{case.name}_slice_xy.svg").write_text(
            dt.slice_product_to_svg(slice_product, width=520),
            encoding="utf-8",
        )
        (output_path / f"{case.name}_tracks_xy.svg").write_text(
            dt.track_product_to_svg(track_product, plane="xy", width=520),
            encoding="utf-8",
        )
        if analytic_slice_xy:
            (output_path / f"{case.name}_analytic_slice_xy.svg").write_text(
                model.quadratic_slice_svg("xy", case.center[2], width=520),
                encoding="utf-8",
            )
        if analytic_projection_xy:
            (output_path / f"{case.name}_analytic_projection_xy.svg").write_text(
                model.quadratic_projection_svg("xy", width=520),
                encoding="utf-8",
            )

    return result


def print_summary(result: dict) -> None:
    print(f"{result['name']}:")
    print(f"  center: {tuple(result['center'])}")
    print(f"  semi_axes: {tuple(round(value, 8) for value in result['semi_axes'])}")
    print(f"  extents: {result['extents']}")
    print(f"  tracks/segments: {result['tracks']} / {result['segments']}")
    print(f"  exact spheroid volume: {result['exact_spheroid_volume']:.8f}")
    print(
        "  tracked spheroid volume: "
        f"{result['tracked_region_volumes'].get(SPHEROID_REGION, 0.0):.8f}"
    )
    print(f"  spheroid volume error: {result['spheroid_volume_error']:.8e}")
    print(f"  total volume error: {result['total_volume_error']:.8e}")
    timings = result["timings"]
    print(
        "  timings [s]: "
        f"build={timings['build_s']:.6f}, "
        f"slice={timings['slice_s']:.6f}, "
        f"track={timings['track_s']:.6f}, "
        f"total={timings['total_s']:.6f}"
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=("all", *CASES.keys()),
        default="centered_sphere",
        help="built-in spheroid case",
    )
    parser.add_argument("--center", nargs=3, type=float, help="custom center")
    parser.add_argument(
        "--axis-ratios", nargs=3, type=float, help="custom spheroid axis ratios"
    )
    parser.add_argument(
        "--directions",
        choices=tuple(DIRECTION_SETS.keys()),
        default="xyz",
        help="tracking direction set",
    )
    parser.add_argument("--spacing", type=float, default=0.05)
    parser.add_argument("--minimum-tracks", type=int, default=1)
    parser.add_argument("--slice-cells", type=int, default=80)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--json", action="store_true", help="emit JSON")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.center is not None or args.axis_ratios is not None:
        center = tuple(args.center or CASES[args.case].center)
        axis_ratios = tuple(args.axis_ratios or CASES[args.case].axis_ratios)
        cases = [SpheroidCase("custom", center=center, axis_ratios=axis_ratios)]
    elif args.case == "all":
        cases = list(CASES.values())
    else:
        cases = [CASES[args.case]]

    results = [
        benchmark_case(
            case,
            directions=DIRECTION_SETS[args.directions],
            spacing=args.spacing,
            minimum_tracks=args.minimum_tracks,
            slice_cells=args.slice_cells,
            repeats=args.repeats,
            output_dir=args.output_dir,
        )
        for case in cases
    ]

    if args.json:
        print(json.dumps(results, indent=2, sort_keys=True))
    else:
        for result in results:
            print_summary(result)


if __name__ == "__main__":
    main()
