"""Hierarchical CSG-style 17x17 PWR assembly model.

This example follows the VERA Core Physics Benchmark progression: pin cell,
2-D lattice, 3-D assembly, spacer grids, and control rods.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


SLICE_RESOLUTION_SCALE = 5
BASE_HOMOGENEOUS_SLICE_CELLS = 32
BASE_PIN_SLICE_CELLS = 80
BASE_ASSEMBLY_SLICE_CELLS = 289
DEFAULT_HOMOGENEOUS_SLICE_CELLS = (
    SLICE_RESOLUTION_SCALE * BASE_HOMOGENEOUS_SLICE_CELLS
)
DEFAULT_PIN_SLICE_CELLS = SLICE_RESOLUTION_SCALE * BASE_PIN_SLICE_CELLS
DEFAULT_ASSEMBLY_SLICE_CELLS = SLICE_RESOLUTION_SCALE * BASE_ASSEMBLY_SLICE_CELLS
DEFAULT_RAY_ANGLES = (0.0, 12.0, 30.0, 45.0, 60.0, 78.0)


OUTSIDE = -1
MODERATOR = 0
FUEL = 1
GAP = 2
CLAD = 3
GUIDE_WATER = 4
GUIDE_TUBE = 5
INSTRUMENT_VOID = 6
INSTRUMENT_THIMBLE = 7
CONTROL_AIC = 8
CONTROL_B4C = 9
CONTROL_CLAD = 10
SPACER_ZIRC = 11
SPACER_INCONEL = 12
BOTTOM_NOZZLE = 13
TOP_NOZZLE = 14


REGION_NAMES: dict[int, str] = {
    MODERATOR: "moderator",
    FUEL: "uo2_fuel",
    GAP: "helium_gap",
    CLAD: "zircaloy_clad",
    GUIDE_WATER: "guide_tube_water",
    GUIDE_TUBE: "zircaloy_guide_tube",
    INSTRUMENT_VOID: "instrument_void",
    INSTRUMENT_THIMBLE: "instrument_thimble",
    CONTROL_AIC: "control_aic",
    CONTROL_B4C: "control_b4c",
    CONTROL_CLAD: "control_clad",
    SPACER_ZIRC: "zircaloy_spacer_grid",
    SPACER_INCONEL: "inconel_spacer_grid",
    BOTTOM_NOZZLE: "bottom_nozzle",
    TOP_NOZZLE: "top_nozzle",
}


REGION_COLORS: dict[int, str] = {
    MODERATOR: "#7fb7d7",
    FUEL: "#d6b34a",
    GAP: "#f3ead7",
    CLAD: "#7a8793",
    GUIDE_WATER: "#9bd0e5",
    GUIDE_TUBE: "#5f6f7a",
    INSTRUMENT_VOID: "#1f2933",
    INSTRUMENT_THIMBLE: "#8f99a3",
    CONTROL_AIC: "#8b5a8c",
    CONTROL_B4C: "#2f4f4f",
    CONTROL_CLAD: "#c0c7cf",
    SPACER_ZIRC: "#426b5f",
    SPACER_INCONEL: "#a2673f",
    BOTTOM_NOZZLE: "#6d7480",
    TOP_NOZZLE: "#4f5965",
}


@dataclass(frozen=True)
class SpacerGrid:
    """Axial spacer-grid homogenization band."""

    name: str
    center: float
    height: float
    region_id: int

    @property
    def bottom(self) -> float:
        return self.center - 0.5 * self.height

    @property
    def top(self) -> float:
        return self.center + 0.5 * self.height

    def contains(self, z: float) -> bool:
        return self.bottom <= z <= self.top


@dataclass(frozen=True)
class PWR17Spec:
    """Centralized VERA-style 17x17 geometry constants in cm."""

    pin_pitch: float = 1.26
    assembly_pitch: float = 21.50
    inter_assembly_half_gap: float = 0.04
    total_height: float = 406.337
    bottom_nozzle_height: float = 6.053
    top_nozzle_height: float = 8.827
    fuel_stack_bottom: float = 11.951
    fuel_stack_height: float = 365.76
    pellet_radius: float = 0.4096
    inner_clad_radius: float = 0.418
    outer_clad_radius: float = 0.475
    inner_guide_tube_radius: float = 0.561
    outer_guide_tube_radius: float = 0.602
    inner_instrument_tube_radius: float = 0.559
    outer_instrument_tube_radius: float = 0.605
    instrument_thimble_inner_radius: float = 0.258
    instrument_thimble_outer_radius: float = 0.382
    control_aic_radius: float = 0.382
    control_b4c_radius: float = 0.373
    control_clad_outer_radius: float = 0.484
    control_aic_height: float = 101.6
    control_b4c_height: float = 259.08
    control_poison_bottom_full_in: float = 17.031
    control_step_size: float = 1.5875
    control_max_steps: int = 230
    spacer_grids: tuple[SpacerGrid, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not self.spacer_grids:
            grids = (
                SpacerGrid("lower_end_grid", 13.884, 3.866, SPACER_INCONEL),
                SpacerGrid("intermediate_grid_1", 75.2, 3.810, SPACER_ZIRC),
                SpacerGrid("intermediate_grid_2", 127.4, 3.810, SPACER_ZIRC),
                SpacerGrid("intermediate_grid_3", 179.6, 3.810, SPACER_ZIRC),
                SpacerGrid("intermediate_grid_4", 231.8, 3.810, SPACER_ZIRC),
                SpacerGrid("intermediate_grid_5", 284.0, 3.810, SPACER_ZIRC),
                SpacerGrid("intermediate_grid_6", 336.2, 3.810, SPACER_ZIRC),
                SpacerGrid("upper_end_grid", 388.2, 3.866, SPACER_INCONEL),
            )
            object.__setattr__(self, "spacer_grids", grids)

    @property
    def number_pins(self) -> int:
        return 17

    @property
    def lattice_width(self) -> float:
        return self.number_pins * self.pin_pitch

    @property
    def fuel_stack_top(self) -> float:
        return self.fuel_stack_bottom + self.fuel_stack_height

    @property
    def fuel_stack_midplane(self) -> float:
        return 0.5 * (self.fuel_stack_bottom + self.fuel_stack_top)

    @property
    def top_nozzle_bottom(self) -> float:
        return self.total_height - self.top_nozzle_height

    @property
    def control_poison_height(self) -> float:
        return self.control_aic_height + self.control_b4c_height

    def pin_center(self, ix: int, iy: int) -> tuple[float, float]:
        x0 = self.inter_assembly_half_gap
        y0 = self.inter_assembly_half_gap
        return (
            x0 + (ix + 0.5) * self.pin_pitch,
            y0 + (iy + 0.5) * self.pin_pitch,
        )


GUIDE_TUBE_LOCATIONS: frozenset[tuple[int, int]] = frozenset(
    {
        (5, 2),
        (8, 2),
        (11, 2),
        (3, 3),
        (13, 3),
        (2, 5),
        (5, 5),
        (8, 5),
        (11, 5),
        (14, 5),
        (2, 8),
        (5, 8),
        (11, 8),
        (14, 8),
        (2, 11),
        (5, 11),
        (8, 11),
        (11, 11),
        (14, 11),
        (3, 13),
        (13, 13),
        (5, 14),
        (8, 14),
        (11, 14),
    }
)

INSTRUMENT_LOCATION = (8, 8)


@dataclass(frozen=True)
class RegionHit:
    """Result of one point query."""

    region_id: int
    name: str
    pin: tuple[int, int] | None = None
    pin_kind: str | None = None


@dataclass(frozen=True)
class SliceResult:
    """Fixed-grid geometry product for inspection and HDF5 export."""

    region_id: np.ndarray
    material_id: np.ndarray
    x_edges: np.ndarray
    y_edges: np.ndarray
    z: float | None
    plane: str = "xy"
    level: str = "assembly3d"

    def to_viz_product(self):
        """Return a pydetran ``SliceProduct`` if the viz helper is available."""

        viz = _load_viz()
        return viz.SliceProduct(
            region_id=self.region_id,
            material_id=self.material_id,
            x_edges=self.x_edges,
            y_edges=self.y_edges,
            plane=self.plane,
        )


@dataclass(frozen=True)
class RayTraceResult:
    """Sampled straight-ray product through an x-y slice.

    This is deliberately shaped like the HDF5 track product: rays have
    entrance/exit points, and each ray owns ordered region-length segments with
    path coordinates.  The segments are sampled, not exact CSG intersections.
    """

    direction: np.ndarray
    track_enter: np.ndarray
    track_exit: np.ndarray
    track_spatial_weight: np.ndarray
    track_length: np.ndarray
    segment_offset: np.ndarray
    segment_count: np.ndarray
    segment_track_index: np.ndarray
    segment_local_index: np.ndarray
    segment_region_id: np.ndarray
    segment_length: np.ndarray
    segment_s_in: np.ndarray
    segment_s_out: np.ndarray
    z: float | None
    angle_degrees: float
    sample_count: int
    plane: str = "xy"
    level: str = "assembly3d"

    @property
    def number_rays(self) -> int:
        return int(self.track_enter.shape[0])

    @property
    def number_segments(self) -> int:
        return int(self.segment_region_id.shape[0])

    def to_viz_track_product(self):
        """Return a pydetran ``TrackProduct`` for existing SVG rendering."""

        viz = _load_viz()
        return viz.TrackProduct(
            directions=viz.DirectionTable(key0=[0], key1=[0], weight=[1.0]),
            tracks=viz.TrackTable(
                direction_index=np.zeros(self.number_rays, dtype=np.int64),
                local_index=np.arange(self.number_rays, dtype=np.int64),
                identifier=np.arange(self.number_rays, dtype=np.int64),
                segment_offset=self.segment_offset,
                segment_count=self.segment_count,
                enter=self.track_enter,
                exit=self.track_exit,
                spatial_weight=self.track_spatial_weight,
                length=self.track_length,
            ),
            segments=viz.SegmentTable(
                track_index=self.segment_track_index,
                local_index=self.segment_local_index,
                region_id=self.segment_region_id,
                length=self.segment_length,
                s_in=self.segment_s_in,
                s_out=self.segment_s_out,
            ),
        )


@dataclass(frozen=True)
class PWR17Model:
    """Hierarchical 17x17 PWR assembly model.

    Levels:
      homogeneous: a single material box
      pin2d: one fuel pin cell in the x-y plane
      assembly2d: a 17x17 x-y slice at ``slice_z``
      assembly3d: full axial classification with spacers and controls
    """

    level: str = "assembly3d"
    control_position: float = 1.0
    include_spacer_grids: bool = True
    homogeneous_region: int = MODERATOR
    spec: PWR17Spec = field(default_factory=PWR17Spec)
    slice_z: float | None = None

    def __post_init__(self) -> None:
        levels = {"homogeneous", "pin2d", "assembly2d", "assembly3d"}
        if self.level not in levels:
            raise ValueError(f"level must be one of {sorted(levels)}")
        if not 0.0 <= self.control_position <= 1.0:
            raise ValueError("control_position must be in [0, 1]")

    @property
    def width_x(self) -> float:
        if self.level == "pin2d":
            return self.spec.pin_pitch
        return self.spec.assembly_pitch

    @property
    def width_y(self) -> float:
        if self.level == "pin2d":
            return self.spec.pin_pitch
        return self.spec.assembly_pitch

    @property
    def height(self) -> float:
        return self.spec.total_height

    @property
    def default_z(self) -> float | None:
        if self.level in {"pin2d", "assembly2d"}:
            return self.slice_z or self.spec.fuel_stack_midplane
        if self.level == "assembly3d":
            return self.slice_z
        return None

    def region_at(self, x: float, y: float, z: float | None = None) -> RegionHit:
        """Classify one point into the deterministic inspection region."""

        z = self.default_z if z is None else z
        if self.level == "homogeneous":
            if self._inside_xy(x, y):
                return self._hit(self.homogeneous_region)
            return self._hit(OUTSIDE)
        if self.level == "pin2d":
            if not (0.0 <= x <= self.spec.pin_pitch and 0.0 <= y <= self.spec.pin_pitch):
                return self._hit(OUTSIDE)
            cx = cy = 0.5 * self.spec.pin_pitch
            return self._classify_fuel_pin(x, y, cx, cy, z)

        if not self._inside_xy(x, y):
            return self._hit(OUTSIDE)
        if self.level == "assembly3d":
            if z is None:
                raise ValueError("assembly3d point queries require z")
            if z < 0.0 or z > self.spec.total_height:
                return self._hit(OUTSIDE)
            if z < self.spec.bottom_nozzle_height:
                return self._hit(BOTTOM_NOZZLE)
            if z >= self.spec.top_nozzle_bottom:
                return self._hit(TOP_NOZZLE)

        pin = self.pin_index(x, y)
        if pin is None:
            return self._maybe_spacer(MODERATOR, z)

        ix, iy = pin
        cx, cy = self.spec.pin_center(ix, iy)
        kind = self.pin_kind(ix, iy)
        if kind == "guide":
            hit = self._classify_guide_tube(x, y, cx, cy, z)
        elif kind == "instrument":
            hit = self._classify_instrument_tube(x, y, cx, cy)
        else:
            hit = self._classify_fuel_pin(x, y, cx, cy, z)
        if hit.region_id == MODERATOR:
            return self._maybe_spacer(MODERATOR, z, pin=pin, pin_kind=kind)
        return hit

    def pin_index(self, x: float, y: float) -> tuple[int, int] | None:
        """Return the 17x17 pin index or ``None`` for inter-assembly gap."""

        x_lattice = x - self.spec.inter_assembly_half_gap
        y_lattice = y - self.spec.inter_assembly_half_gap
        if x_lattice < 0.0 or y_lattice < 0.0:
            return None
        if x_lattice >= self.spec.lattice_width or y_lattice >= self.spec.lattice_width:
            return None
        ix = min(int(x_lattice / self.spec.pin_pitch), self.spec.number_pins - 1)
        iy = min(int(y_lattice / self.spec.pin_pitch), self.spec.number_pins - 1)
        return ix, iy

    def pin_kind(self, ix: int, iy: int) -> str:
        if (ix, iy) == INSTRUMENT_LOCATION:
            return "instrument"
        if (ix, iy) in GUIDE_TUBE_LOCATIONS:
            return "guide"
        return "fuel"

    def control_poison_bounds(self) -> tuple[float, float]:
        """Return bottom/top of the poison stack for the insertion fraction."""

        withdrawn_steps = (1.0 - self.control_position) * self.spec.control_max_steps
        bottom = (
            self.spec.control_poison_bottom_full_in
            + withdrawn_steps * self.spec.control_step_size
        )
        return bottom, bottom + self.spec.control_poison_height

    def spacer_grid_at(self, z: float | None) -> SpacerGrid | None:
        if z is None or not self.include_spacer_grids:
            return None
        for grid in self.spec.spacer_grids:
            if grid.contains(z):
                return grid
        return None

    def sample_slice(self, nx: int, ny: int, z: float | None = None) -> SliceResult:
        """Gridify the analytic model into a reusable slice product."""

        if nx <= 0 or ny <= 0:
            raise ValueError("nx and ny must be positive")
        z = self.default_z if z is None else z
        x_edges = np.linspace(0.0, self.width_x, nx + 1)
        y_edges = np.linspace(0.0, self.width_y, ny + 1)
        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
        region_id = np.empty((ny, nx), dtype=np.int64)
        for j, y in enumerate(y_centers):
            for i, x in enumerate(x_centers):
                region_id[j, i] = self.region_at(float(x), float(y), z).region_id
        return SliceResult(
            region_id=region_id,
            material_id=region_id.copy(),
            x_edges=x_edges,
            y_edges=y_edges,
            z=z,
            level=self.level,
        )

    def trace_rays(
        self,
        z: float | None = None,
        *,
        angle_degrees: float = 0.0,
        number_rays: int = 57,
        sample_count: int = 1200,
    ) -> RayTraceResult:
        """Sample straight rays across an x-y slice and merge region spans."""

        if number_rays <= 0:
            raise ValueError("number_rays must be positive")
        if sample_count <= 0:
            raise ValueError("sample_count must be positive")
        z = self.default_z if z is None else z
        if z is None and self.level == "assembly3d":
            z = self.spec.fuel_stack_midplane

        theta = np.deg2rad(angle_degrees)
        direction = np.array([np.cos(theta), np.sin(theta)], dtype=float)
        normal = np.array([-direction[1], direction[0]], dtype=float)
        corners = np.array(
            [
                [0.0, 0.0],
                [self.width_x, 0.0],
                [0.0, self.width_y],
                [self.width_x, self.width_y],
            ],
            dtype=float,
        )
        projection = corners @ normal
        offset_edges = np.linspace(
            float(np.min(projection)), float(np.max(projection)), number_rays + 1
        )
        offsets = 0.5 * (offset_edges[:-1] + offset_edges[1:])
        ray_width = float(offset_edges[1] - offset_edges[0])

        track_enter: list[list[float]] = []
        track_exit: list[list[float]] = []
        track_spatial_weight: list[float] = []
        track_length: list[float] = []
        segment_offset: list[int] = []
        segment_count: list[int] = []
        segment_track_index: list[int] = []
        segment_local_index: list[int] = []
        segment_region_id: list[int] = []
        segment_length: list[float] = []
        segment_s_in: list[float] = []
        segment_s_out: list[float] = []

        z_coord = 0.0 if z is None else float(z)
        for offset in offsets:
            clipped = _clip_line_to_box(
                origin=normal * offset,
                direction=direction,
                width=self.width_x,
                height=self.width_y,
            )
            if clipped is None:
                continue
            p0, p1 = clipped
            length = float(np.linalg.norm(p1 - p0))
            if length <= 0.0:
                continue

            ray_index = len(track_enter)
            track_enter.append([float(p0[0]), float(p0[1]), z_coord])
            track_exit.append([float(p1[0]), float(p1[1]), z_coord])
            track_spatial_weight.append(ray_width)
            track_length.append(length)
            segment_offset.append(len(segment_region_id))

            edges = np.linspace(0.0, length, sample_count + 1)
            centers = 0.5 * (edges[:-1] + edges[1:])
            regions = []
            for s in centers:
                point = p0 + direction * s
                regions.append(self.region_at(float(point[0]), float(point[1]), z).region_id)

            start = 0
            local_index = 0
            for i in range(1, len(regions) + 1):
                if i < len(regions) and regions[i] == regions[start]:
                    continue
                s_in = float(edges[start])
                s_out = float(edges[i])
                segment_track_index.append(ray_index)
                segment_local_index.append(local_index)
                segment_region_id.append(int(regions[start]))
                segment_length.append(s_out - s_in)
                segment_s_in.append(s_in)
                segment_s_out.append(s_out)
                start = i
                local_index += 1
            segment_count.append(local_index)

        return RayTraceResult(
            direction=np.array([direction[0], direction[1], 0.0], dtype=float),
            track_enter=np.asarray(track_enter, dtype=float),
            track_exit=np.asarray(track_exit, dtype=float),
            track_spatial_weight=np.asarray(track_spatial_weight, dtype=float),
            track_length=np.asarray(track_length, dtype=float),
            segment_offset=np.asarray(segment_offset, dtype=np.int64),
            segment_count=np.asarray(segment_count, dtype=np.int64),
            segment_track_index=np.asarray(segment_track_index, dtype=np.int64),
            segment_local_index=np.asarray(segment_local_index, dtype=np.int64),
            segment_region_id=np.asarray(segment_region_id, dtype=np.int64),
            segment_length=np.asarray(segment_length, dtype=float),
            segment_s_in=np.asarray(segment_s_in, dtype=float),
            segment_s_out=np.asarray(segment_s_out, dtype=float),
            z=z,
            angle_degrees=angle_degrees,
            sample_count=sample_count,
            level=self.level,
        )

    def _inside_xy(self, x: float, y: float) -> bool:
        return 0.0 <= x <= self.width_x and 0.0 <= y <= self.width_y

    def _active_fuel_z(self, z: float | None) -> bool:
        if z is None:
            return True
        return self.spec.fuel_stack_bottom <= z <= self.spec.fuel_stack_top

    def _radius(self, x: float, y: float, cx: float, cy: float) -> float:
        return float(np.hypot(x - cx, y - cy))

    def _classify_fuel_pin(
        self, x: float, y: float, cx: float, cy: float, z: float | None
    ) -> RegionHit:
        r = self._radius(x, y, cx, cy)
        if self._active_fuel_z(z):
            if r <= self.spec.pellet_radius:
                return self._hit(FUEL)
            if r <= self.spec.inner_clad_radius:
                return self._hit(GAP)
            if r <= self.spec.outer_clad_radius:
                return self._hit(CLAD)
            return self._hit(MODERATOR)
        if r <= self.spec.outer_clad_radius:
            return self._hit(CLAD)
        return self._hit(MODERATOR)

    def _classify_guide_tube(
        self, x: float, y: float, cx: float, cy: float, z: float | None
    ) -> RegionHit:
        r = self._radius(x, y, cx, cy)
        if z is not None:
            bottom, top = self.control_poison_bounds()
            if bottom <= z <= top:
                if z <= bottom + self.spec.control_aic_height:
                    poison_radius = self.spec.control_aic_radius
                    poison_region = CONTROL_AIC
                else:
                    poison_radius = self.spec.control_b4c_radius
                    poison_region = CONTROL_B4C
                if r <= poison_radius:
                    return self._hit(poison_region)
                if r <= self.spec.control_clad_outer_radius:
                    return self._hit(CONTROL_CLAD)
        if r <= self.spec.inner_guide_tube_radius:
            return self._hit(GUIDE_WATER)
        if r <= self.spec.outer_guide_tube_radius:
            return self._hit(GUIDE_TUBE)
        return self._hit(MODERATOR)

    def _classify_instrument_tube(
        self, x: float, y: float, cx: float, cy: float
    ) -> RegionHit:
        r = self._radius(x, y, cx, cy)
        if r <= self.spec.instrument_thimble_inner_radius:
            return self._hit(INSTRUMENT_VOID)
        if r <= self.spec.instrument_thimble_outer_radius:
            return self._hit(INSTRUMENT_THIMBLE)
        if r <= self.spec.inner_instrument_tube_radius:
            return self._hit(GUIDE_WATER)
        if r <= self.spec.outer_instrument_tube_radius:
            return self._hit(GUIDE_TUBE)
        return self._hit(MODERATOR)

    def _maybe_spacer(
        self,
        base_region: int,
        z: float | None,
        pin: tuple[int, int] | None = None,
        pin_kind: str | None = None,
    ) -> RegionHit:
        grid = self.spacer_grid_at(z)
        if grid is not None and base_region == MODERATOR:
            return self._hit(grid.region_id, pin=pin, pin_kind=pin_kind)
        return self._hit(base_region, pin=pin, pin_kind=pin_kind)

    def _hit(
        self,
        region_id: int,
        pin: tuple[int, int] | None = None,
        pin_kind: str | None = None,
    ) -> RegionHit:
        if region_id == OUTSIDE:
            return RegionHit(region_id=OUTSIDE, name="outside")
        return RegionHit(
            region_id=region_id,
            name=REGION_NAMES[region_id],
            pin=pin,
            pin_kind=pin_kind,
        )


def homogeneous_model(region_id: int = MODERATOR) -> PWR17Model:
    return PWR17Model(level="homogeneous", homogeneous_region=region_id)


def fuel_pin_2d_model() -> PWR17Model:
    return PWR17Model(level="pin2d")


def assembly_2d_model(control_position: float = 1.0) -> PWR17Model:
    return PWR17Model(level="assembly2d", control_position=control_position)


def assembly_3d_model(
    control_position: float = 1.0, include_spacer_grids: bool = True
) -> PWR17Model:
    return PWR17Model(
        level="assembly3d",
        control_position=control_position,
        include_spacer_grids=include_spacer_grids,
    )


def slice_to_svg(
    product: SliceResult,
    *,
    width: int = 900,
    height: int | None = None,
    title: str | None = None,
) -> str:
    """Render a sampled slice with pydetran's backend-neutral SVG helper."""

    viz = _load_viz()
    title = title or f"PWR17 {product.level} slice"
    return viz.slice_product_to_svg(
        product.to_viz_product(),
        width=width,
        height=height,
        region_colors=REGION_COLORS,
        title=title,
    )


def slice_to_png(
    product: SliceResult,
    filename: str | Path,
    *,
    region_colors: Mapping[int, str] | None = None,
) -> None:
    """Write a fixed-grid slice as a compact raster image."""

    try:
        from PIL import Image
    except ImportError as exc:
        raise ImportError("slice PNG output requires Pillow") from exc

    colors = region_colors or REGION_COLORS
    rgb = np.empty(product.region_id.shape + (3,), dtype=np.uint8)
    for region_id in np.unique(product.region_id):
        rgb[product.region_id == region_id] = _rgb_for_region(int(region_id), colors)
    image = Image.fromarray(np.flipud(rgb), mode="RGB")
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def ray_trace_to_svg(
    product: RayTraceResult,
    *,
    width: int = 900,
    height: int | None = None,
    stroke_width: float = 1.5,
    title: str | None = None,
) -> str:
    """Render sampled rays as region-colored vector segments."""

    viz = _load_viz()
    title = title or f"PWR17 {product.level} ray trace"
    return viz.track_product_to_svg(
        product.to_viz_track_product(),
        width=width,
        height=height,
        stroke_width=stroke_width,
        region_colors=REGION_COLORS,
        title=title,
    )


def write_hdf5(
    filename: str | Path,
    model: PWR17Model,
    slices: Mapping[str, SliceResult],
    rays: Mapping[str, RayTraceResult] | None = None,
) -> None:
    """Write source metadata plus sampled slice/ray products to HDF5."""

    try:
        import h5py
    except ImportError as exc:
        raise ImportError("write_hdf5 requires h5py") from exc

    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as h5:
        geometry = h5.create_group("geometry")
        source = geometry.create_group("source")
        source.attrs["backend"] = "python-csg"
        source.attrs["model"] = "pwr17_assembly"
        source.attrs["level"] = model.level
        source.attrs["schema_version"] = 1
        source.attrs["control_position"] = model.control_position
        source.attrs["include_spacer_grids"] = int(model.include_spacer_grids)
        source.attrs["benchmark_family"] = (
            "VERA Core Physics Benchmark Progression Problems"
        )
        source.attrs["benchmark_note"] = (
            "VERA-style scaffold for inspection, not a certified benchmark input"
        )

        constants = source.create_group("constants")
        for key, value in model.spec.__dict__.items():
            if key == "spacer_grids":
                continue
            constants.attrs[key] = value
        grids = source.create_group("spacer_grids")
        grids.create_dataset(
            "name", data=[g.name for g in model.spec.spacer_grids], dtype=string_dtype
        )
        grids.create_dataset("center", data=[g.center for g in model.spec.spacer_grids])
        grids.create_dataset("height", data=[g.height for g in model.spec.spacer_grids])
        grids.create_dataset(
            "region_id", data=[g.region_id for g in model.spec.spacer_grids]
        )

        regions = geometry.create_group("regions")
        ids = np.array(sorted(REGION_NAMES), dtype=np.int64)
        regions.create_dataset("id", data=ids)
        regions.create_dataset(
            "name", data=[REGION_NAMES[int(i)] for i in ids], dtype=string_dtype
        )
        regions.create_dataset("material_id", data=ids)
        regions.create_dataset(
            "material_name",
            data=[REGION_NAMES[int(i)] for i in ids],
            dtype=string_dtype,
        )

        products = geometry.create_group("products")
        slices_root = products.create_group("slices")
        slices_root.attrs["schema_version"] = 1
        for name, product in slices.items():
            group = slices_root.create_group(name)
            group.attrs["plane"] = product.plane
            group.attrs["level"] = product.level
            if product.z is not None:
                group.attrs["z"] = product.z
            group.create_dataset("region_id", data=product.region_id)
            group.create_dataset("material_id", data=product.material_id)
            group.create_dataset("x_edges", data=product.x_edges)
            group.create_dataset("y_edges", data=product.y_edges)

        if rays:
            rays_root = products.create_group("rays")
            rays_root.attrs["schema_version"] = 1
            for name, product in rays.items():
                group = rays_root.create_group(name)
                group.attrs["plane"] = product.plane
                group.attrs["level"] = product.level
                group.attrs["angle_degrees"] = product.angle_degrees
                group.attrs["sample_count"] = product.sample_count
                if product.z is not None:
                    group.attrs["z"] = product.z
                group.create_dataset("direction", data=product.direction)
                tracks = group.create_group("tracks")
                tracks.create_dataset("enter", data=product.track_enter)
                tracks.create_dataset("exit", data=product.track_exit)
                tracks.create_dataset("spatial_weight", data=product.track_spatial_weight)
                tracks.create_dataset("length", data=product.track_length)
                tracks.create_dataset("segment_offset", data=product.segment_offset)
                tracks.create_dataset("segment_count", data=product.segment_count)
                segments = group.create_group("segments")
                segments.create_dataset("track_index", data=product.segment_track_index)
                segments.create_dataset("local_index", data=product.segment_local_index)
                segments.create_dataset("region_id", data=product.segment_region_id)
                segments.create_dataset("length", data=product.segment_length)
                segments.create_dataset("s_in", data=product.segment_s_in)
                segments.create_dataset("s_out", data=product.segment_s_out)


def example_slices(
    control_position: float = 1.0,
    nx: int = DEFAULT_ASSEMBLY_SLICE_CELLS,
    ny: int = DEFAULT_ASSEMBLY_SLICE_CELLS,
) -> dict[str, SliceResult]:
    """Build a small hierarchy of slices from simple to detailed."""

    spec = PWR17Spec()
    slices: dict[str, SliceResult] = {}
    slices["homogeneous"] = homogeneous_model().sample_slice(
        DEFAULT_HOMOGENEOUS_SLICE_CELLS,
        DEFAULT_HOMOGENEOUS_SLICE_CELLS,
    )
    slices["fuel_pin_2d"] = fuel_pin_2d_model().sample_slice(
        DEFAULT_PIN_SLICE_CELLS,
        DEFAULT_PIN_SLICE_CELLS,
    )
    slices["assembly_2d"] = assembly_2d_model(control_position).sample_slice(nx, ny)
    assembly = assembly_3d_model(control_position)
    slices["assembly_3d_midplane"] = assembly.sample_slice(nx, ny, spec.fuel_stack_midplane)
    slices["assembly_3d_spacer_grid"] = assembly.sample_slice(nx, ny, 75.2)
    return slices


def example_ray_traces(
    control_position: float = 1.0,
    *,
    number_rays: int = 57,
    sample_count: int = 1200,
    angle_degrees: float | Sequence[float] = DEFAULT_RAY_ANGLES,
) -> dict[str, RayTraceResult]:
    """Build ray products that expose track-like CSG inspection data."""

    spec = PWR17Spec()
    assembly = assembly_3d_model(control_position)
    ray_products: dict[str, RayTraceResult] = {}
    planes = (
        ("assembly_3d_midplane_rays", spec.fuel_stack_midplane),
        ("assembly_3d_spacer_grid_rays", 75.2),
    )
    for angle in _angle_sequence(angle_degrees):
        token = _angle_token(angle)
        for name, z in planes:
            ray_products[f"{name}_{token}"] = assembly.trace_rays(
                z,
                angle_degrees=angle,
                number_rays=number_rays,
                sample_count=sample_count,
            )
    return ray_products


def write_example_outputs(
    output_dir: str | Path,
    *,
    control_position: float = 1.0,
    nx: int = DEFAULT_ASSEMBLY_SLICE_CELLS,
    ny: int = DEFAULT_ASSEMBLY_SLICE_CELLS,
    rays: bool = True,
    number_rays: int = 57,
    ray_samples: int = 1200,
    ray_angle: float | Sequence[float] = DEFAULT_RAY_ANGLES,
    hdf5: bool = True,
    slice_format: str = "png",
) -> list[Path]:
    """Write slice images, ray SVGs, and optional HDF5 outputs."""

    if slice_format not in {"png", "svg", "both"}:
        raise ValueError("slice_format must be 'png', 'svg', or 'both'")

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    slices = example_slices(control_position, nx, ny)
    ray_products = (
        example_ray_traces(
            control_position,
            number_rays=number_rays,
            sample_count=ray_samples,
            angle_degrees=ray_angle,
        )
        if rays
        else {}
    )
    written: list[Path] = []
    for name, product in slices.items():
        if slice_format in {"png", "both"}:
            path = output / f"{name}.png"
            slice_to_png(product, path)
            written.append(path)
        if slice_format in {"svg", "both"}:
            path = output / f"{name}.svg"
            path.write_text(slice_to_svg(product), encoding="utf-8")
            written.append(path)
    for name, product in ray_products.items():
        path = output / f"{name}.svg"
        path.write_text(ray_trace_to_svg(product), encoding="utf-8")
        written.append(path)
    if hdf5:
        h5_path = output / "pwr17_csg_products.h5"
        write_hdf5(h5_path, assembly_3d_model(control_position), slices, ray_products)
        written.append(h5_path)
    return written


def _clip_line_to_box(
    origin: np.ndarray, direction: np.ndarray, width: float, height: float
) -> tuple[np.ndarray, np.ndarray] | None:
    t_min = -np.inf
    t_max = np.inf
    for coordinate, slope, low, high in (
        (origin[0], direction[0], 0.0, width),
        (origin[1], direction[1], 0.0, height),
    ):
        if abs(float(slope)) < 1.0e-14:
            if coordinate < low or coordinate > high:
                return None
            continue
        t0 = (low - coordinate) / slope
        t1 = (high - coordinate) / slope
        t_low = min(t0, t1)
        t_high = max(t0, t1)
        t_min = max(t_min, t_low)
        t_max = min(t_max, t_high)
    if not np.isfinite(t_min) or not np.isfinite(t_max) or t_max <= t_min:
        return None
    return origin + t_min * direction, origin + t_max * direction


def _angle_sequence(angle_degrees: float | Sequence[float]) -> tuple[float, ...]:
    if isinstance(angle_degrees, (str, bytes)):
        raise TypeError("angle_degrees must be a float or sequence of floats")
    try:
        return (float(angle_degrees),)
    except TypeError:
        return tuple(float(angle) for angle in angle_degrees)


def _angle_token(angle_degrees: float) -> str:
    text = f"{angle_degrees:.3f}".rstrip("0").rstrip(".")
    text = text.replace("-", "m").replace(".", "p")
    return f"angle_{text}"


def _rgb_for_region(region_id: int, colors: Mapping[int, str]) -> tuple[int, int, int]:
    color = colors.get(region_id, "#ffffff" if region_id == OUTSIDE else "#000000")
    return _hex_color_to_rgb(color)


def _hex_color_to_rgb(color: str) -> tuple[int, int, int]:
    color = color.strip()
    if color.startswith("#"):
        color = color[1:]
    if len(color) != 6:
        raise ValueError(f"expected #rrggbb color, got {color!r}")
    return int(color[0:2], 16), int(color[2:4], 16), int(color[4:6], 16)


def _load_viz():
    try:
        import detran_viz

        return detran_viz
    except ImportError:
        try:
            from detran import detran_viz

            return detran_viz
        except ImportError:
            pass
    try:
        import detran as dt

        return dt
    except ImportError as exc:
        raise ImportError(
            "CSG SVG output requires pydetran's detran_viz helpers. "
            "Set PYTHONPATH to the pydetran source/build or install pydetran."
        ) from exc


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build VERA-style PWR17 CSG inspection products."
    )
    parser.add_argument(
        "--output-dir",
        default="pwr17_outputs",
        help="directory for SVG and HDF5 products",
    )
    parser.add_argument(
        "--control-position",
        type=float,
        default=1.0,
        help="control insertion fraction, 0 withdrawn and 1 fully inserted",
    )
    parser.add_argument(
        "--nx",
        type=int,
        default=DEFAULT_ASSEMBLY_SLICE_CELLS,
        help="slice cells in x",
    )
    parser.add_argument(
        "--ny",
        type=int,
        default=DEFAULT_ASSEMBLY_SLICE_CELLS,
        help="slice cells in y",
    )
    parser.add_argument(
        "--no-rays",
        action="store_true",
        help="skip sampled ray SVG/HDF5 products",
    )
    parser.add_argument("--number-rays", type=int, default=57, help="inspection rays")
    parser.add_argument(
        "--ray-samples",
        type=int,
        default=1200,
        help="point samples per inspection ray",
    )
    parser.add_argument(
        "--ray-angle",
        "--ray-angles",
        nargs="+",
        type=float,
        default=list(DEFAULT_RAY_ANGLES),
        help="ray angle(s) in degrees in the x-y plane",
    )
    parser.add_argument(
        "--slice-format",
        choices=("png", "svg", "both"),
        default="png",
        help="slice image format; png is best for high-resolution pixel maps",
    )
    parser.add_argument(
        "--no-hdf5",
        action="store_true",
        help="write SVG files only, without the optional HDF5 product",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    paths = write_example_outputs(
        args.output_dir,
        control_position=args.control_position,
        nx=args.nx,
        ny=args.ny,
        rays=not args.no_rays,
        number_rays=args.number_rays,
        ray_samples=args.ray_samples,
        ray_angle=args.ray_angle,
        hdf5=not args.no_hdf5,
        slice_format=args.slice_format,
    )
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
