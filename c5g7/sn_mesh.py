"""Automated Cartesian pin meshes for C5G7 SN/DD studies.

The first concrete target is the 7x7, volume-preserving circular-pin
approximation used by the legacy C5G7 SN drivers.  The generator below keeps
the same D4-symmetric 13-cell central-material mask but computes the widths
from a one-dimensional objective instead of baking the historical constants
into the mesh construction path.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable, Sequence


PIN_PITCH = 1.26
PIN_RADIUS = 0.54
ASSEMBLY_DIMENSION = 17
CORE_DIMENSION = 3

PIN_NAMES = {
  0: "UO2 fuel",
  1: "MOX 4.3% fuel",
  2: "MOX 7.0% fuel",
  3: "MOX 8.7% fuel",
  4: "guide tube",
  5: "fission chamber",
  6: "moderator",
}

CENTRAL_MATERIAL_MASK_7 = (
  0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0,
  0, 0, 1, 1, 1, 0, 0,
  0, 1, 1, 1, 1, 1, 0,
  0, 0, 1, 1, 1, 0, 0,
  0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0,
)


@dataclass(frozen=True)
class PinMeshRecipe:
  """A reusable C5G7 pin mesh recipe."""

  name: str
  number: int
  pitch: float
  radius: float
  widths: tuple[float, ...]
  edges: tuple[float, ...]
  central_material_mask: tuple[int, ...]
  central_cell_count: int
  central_area_target: float
  central_area_mesh: float
  central_area_error: float
  moderator_area_target: float
  moderator_area_mesh: float
  moderator_area_error: float
  central_overlap_area: float
  false_central_area: float
  missed_central_area: float
  symmetric_difference_area: float
  symmetric_difference_fraction: float
  optimized_l_over_r: float
  delta_over_r: float
  legacy_l_over_r: float | None = None
  legacy_symmetric_difference_fraction: float | None = None

  def to_dict(self) -> dict:
    return asdict(self)


@dataclass(frozen=True)
class MaterialAreaCorrection:
  """Area inventory correction for one material in a center-sampled mesh."""

  material: int
  material_name: str
  target_area: float
  assigned_area: float
  target_fraction: float
  assigned_fraction: float
  area_error: float
  fraction_error: float
  density_correction: float | None
  conservative: bool


@dataclass(frozen=True)
class CenterSampledMeshRecipe:
  """Problem-specific center-sampled C5G7 pin-mesh recipe."""

  name: str
  number: int
  pitch: float
  radius: float
  widths: tuple[float, ...]
  edges: tuple[float, ...]
  pins_per_side: int
  pin_count: int
  total_area: float
  correction_scope: str
  material_area_corrections: tuple[MaterialAreaCorrection, ...]

  @property
  def conservative(self) -> bool:
    return all(entry.conservative for entry in self.material_area_corrections)

  def density_corrections(self) -> dict[int, float]:
    """Return finite density corrections keyed by material id."""

    return {
      entry.material: float(entry.density_correction)
      for entry in self.material_area_corrections
      if entry.density_correction is not None
    }

  def to_dict(self) -> dict:
    entries = [asdict(entry) for entry in self.material_area_corrections]
    return {
      "name": self.name,
      "number": self.number,
      "pitch": self.pitch,
      "radius": self.radius,
      "widths": self.widths,
      "edges": self.edges,
      "pins_per_side": self.pins_per_side,
      "pin_count": self.pin_count,
      "total_area": self.total_area,
      "correction_scope": self.correction_scope,
      "center_sampled": True,
      "density_corrected": True,
      "conservative_density_correction": self.conservative,
      "material_area_corrections": entries,
      "density_correction_by_material": {
        str(entry.material): entry.density_correction
        for entry in self.material_area_corrections
        if entry.density_correction is not None
      },
      "area_assignment": (
        "Each Cartesian subcell inherits the material at its center. "
        "Macroscopic cross sections are scaled per material so assigned mesh "
        "area times corrected density matches the known analytic material area."
      ),
    }


def _legendre_nodes_weights(n: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
  """Return n-point Gauss-Legendre nodes and weights on [-1, 1]."""

  nodes = [0.0] * n
  weights = [0.0] * n
  m = (n + 1) // 2
  for i in range(m):
    z = math.cos(math.pi * (i + 0.75) / (n + 0.5))
    while True:
      p1 = 1.0
      p2 = 0.0
      for j in range(1, n + 1):
        p3 = p2
        p2 = p1
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
      pp = n * (z * p1 - p2) / (z * z - 1.0)
      z_next = z - p1 / pp
      if abs(z_next - z) < 1.0e-15:
        z = z_next
        break
      z = z_next
    nodes[i] = -z
    nodes[n - 1 - i] = z
    weight = 2.0 / ((1.0 - z * z) * pp * pp)
    weights[i] = weight
    weights[n - 1 - i] = weight
  return tuple(nodes), tuple(weights)


_GL_NODES, _GL_WEIGHTS = _legendre_nodes_weights(24)


def _edges_from_widths(widths: Sequence[float]) -> tuple[float, ...]:
  edges = [0.0]
  for width in widths:
    edges.append(edges[-1] + float(width))
  return tuple(edges)


def _circle_rectangle_overlap(
  x0: float,
  x1: float,
  y0: float,
  y1: float,
  *,
  pitch: float,
  radius: float,
) -> float:
  """Area of the centered circle inside one axis-aligned rectangle."""

  center = 0.5 * pitch
  a = max(x0, center - radius)
  b = min(x1, center + radius)
  if b <= a or y1 <= center - radius or y0 >= center + radius:
    return 0.0

  mid = 0.5 * (a + b)
  half = 0.5 * (b - a)
  integral = 0.0
  for node, weight in zip(_GL_NODES, _GL_WEIGHTS):
    x = mid + half * node
    radicand = radius * radius - (x - center) * (x - center)
    if radicand <= 0.0:
      continue
    height = math.sqrt(radicand)
    lower = max(y0, center - height)
    upper = min(y1, center + height)
    if upper > lower:
      integral += weight * (upper - lower)
  return half * integral


def _central_stats(
  widths: Sequence[float],
  mask: Sequence[int],
  *,
  pitch: float,
  radius: float,
) -> dict:
  number = len(widths)
  edges = _edges_from_widths(widths)
  central_area = 0.0
  overlap = 0.0
  for j in range(number):
    dy = edges[j + 1] - edges[j]
    for i in range(number):
      if not mask[j * number + i]:
        continue
      dx = edges[i + 1] - edges[i]
      central_area += dx * dy
      overlap += _circle_rectangle_overlap(
        edges[i],
        edges[i + 1],
        edges[j],
        edges[j + 1],
        pitch=pitch,
        radius=radius,
      )
  target = math.pi * radius * radius
  false_central = central_area - overlap
  missed_central = target - overlap
  symmetric_difference = false_central + missed_central
  return {
    "central_area_mesh": central_area,
    "central_overlap_area": overlap,
    "false_central_area": false_central,
    "missed_central_area": missed_central,
    "symmetric_difference_area": symmetric_difference,
    "symmetric_difference_fraction": symmetric_difference / target,
  }


def _one_step_widths(l_over_r: float, *, pitch: float, radius: float) -> tuple[float, ...]:
  """Widths for the 7x7 one-step, D4-symmetric circular-pin surrogate."""

  length = l_over_r * radius
  delta = (math.pi - l_over_r * l_over_r) * radius / (2.0 * l_over_r)
  moderator = 0.5 * pitch - 0.5 * length - delta
  quarter = 0.25 * length
  half = 0.5 * length
  return (moderator, delta, quarter, half, quarter, delta, moderator)


def _objective(l_over_r: float, *, pitch: float, radius: float) -> float:
  widths = _one_step_widths(l_over_r, pitch=pitch, radius=radius)
  if min(widths) <= 0.0:
    return float("inf")
  return _central_stats(
    widths,
    CENTRAL_MATERIAL_MASK_7,
    pitch=pitch,
    radius=radius,
  )["symmetric_difference_fraction"]


def _optimize_l_over_r(*, pitch: float, radius: float) -> float:
  """Golden-section minimization for the area-preserving 7x7 shape."""

  lower = 1.0
  upper = math.sqrt(math.pi) - 1.0e-12
  while _one_step_widths(upper, pitch=pitch, radius=radius)[0] <= 0.0:
    upper -= 1.0e-3
  phi = 0.5 * (1.0 + math.sqrt(5.0))
  c = upper - (upper - lower) / phi
  d = lower + (upper - lower) / phi
  fc = _objective(c, pitch=pitch, radius=radius)
  fd = _objective(d, pitch=pitch, radius=radius)
  for _ in range(96):
    if fc < fd:
      upper = d
      d = c
      fd = fc
      c = upper - (upper - lower) / phi
      fc = _objective(c, pitch=pitch, radius=radius)
    else:
      lower = c
      c = d
      fc = fd
      d = lower + (upper - lower) / phi
      fd = _objective(d, pitch=pitch, radius=radius)
  return 0.5 * (lower + upper)


def optimized_pin_recipe_7(
  *,
  pitch: float = PIN_PITCH,
  radius: float = PIN_RADIUS,
) -> PinMeshRecipe:
  """Generate the automated 7x7 C5G7 pin recipe."""

  l_over_r = _optimize_l_over_r(pitch=pitch, radius=radius)
  widths = _one_step_widths(l_over_r, pitch=pitch, radius=radius)
  edges = _edges_from_widths(widths)
  stats = _central_stats(
    widths,
    CENTRAL_MATERIAL_MASK_7,
    pitch=pitch,
    radius=radius,
  )
  target = math.pi * radius * radius
  total = pitch * pitch
  legacy_l_over_r = 1.60496875
  legacy_stats = _central_stats(
    _one_step_widths(legacy_l_over_r, pitch=pitch, radius=radius),
    CENTRAL_MATERIAL_MASK_7,
    pitch=pitch,
    radius=radius,
  )
  return PinMeshRecipe(
    name="optimized_7x7_one_step",
    number=7,
    pitch=pitch,
    radius=radius,
    widths=tuple(widths),
    edges=tuple(edges),
    central_material_mask=CENTRAL_MATERIAL_MASK_7,
    central_cell_count=sum(CENTRAL_MATERIAL_MASK_7),
    central_area_target=target,
    central_area_mesh=stats["central_area_mesh"],
    central_area_error=stats["central_area_mesh"] - target,
    moderator_area_target=total - target,
    moderator_area_mesh=total - stats["central_area_mesh"],
    moderator_area_error=(total - stats["central_area_mesh"]) - (total - target),
    central_overlap_area=stats["central_overlap_area"],
    false_central_area=stats["false_central_area"],
    missed_central_area=stats["missed_central_area"],
    symmetric_difference_area=stats["symmetric_difference_area"],
    symmetric_difference_fraction=stats["symmetric_difference_fraction"],
    optimized_l_over_r=l_over_r,
    delta_over_r=widths[1] / radius,
    legacy_l_over_r=legacy_l_over_r,
    legacy_symmetric_difference_fraction=legacy_stats["symmetric_difference_fraction"],
  )


def uniform_pin_widths(number: int, *, pitch: float = PIN_PITCH) -> tuple[float, ...]:
  """Return equal-width pin-cell intervals for a center-sampled mesh."""

  number = int(number)
  if number < 1:
    raise ValueError("number must be positive")
  width = float(pitch) / number
  return tuple(width for _ in range(number))


def _pin_true_material_areas(
  central_material: int,
  *,
  pitch: float,
  radius: float,
  moderator_material: int = 6,
) -> dict[int, float]:
  total = pitch * pitch
  if central_material == moderator_material:
    return {moderator_material: total}
  central_area = math.pi * radius * radius
  return {
    int(central_material): central_area,
    moderator_material: total - central_area,
  }


def center_sampled_material_map_for_pin(
  central_material: int,
  *,
  widths: Sequence[float] | None = None,
  number: int | None = None,
  pitch: float = PIN_PITCH,
  radius: float = PIN_RADIUS,
  moderator_material: int = 6,
) -> list[int]:
  """Return a row-major pin material map assigned by subcell centers."""

  if widths is None:
    if number is None:
      number = 7
    widths = uniform_pin_widths(number, pitch=pitch)
  edges = _edges_from_widths(widths)
  center = 0.5 * pitch
  materials = []
  for j in range(len(widths)):
    y = 0.5 * (edges[j] + edges[j + 1])
    for i in range(len(widths)):
      x = 0.5 * (edges[i] + edges[i + 1])
      if (
        central_material != moderator_material
        and (x - center) * (x - center) + (y - center) * (y - center) <= radius * radius
      ):
        materials.append(int(central_material))
      else:
        materials.append(moderator_material)
  return materials


def _pin_assigned_material_areas(
  central_material: int,
  *,
  widths: Sequence[float],
  pitch: float,
  radius: float,
  moderator_material: int = 6,
) -> dict[int, float]:
  edges = _edges_from_widths(widths)
  materials = center_sampled_material_map_for_pin(
    central_material,
    widths=widths,
    pitch=pitch,
    radius=radius,
    moderator_material=moderator_material,
  )
  areas: dict[int, float] = {}
  number = len(widths)
  for j in range(number):
    dy = edges[j + 1] - edges[j]
    for i in range(number):
      dx = edges[i + 1] - edges[i]
      material = int(materials[j * number + i])
      areas[material] = areas.get(material, 0.0) + dx * dy
  return areas


def _pin_counts(pin_map: Sequence[int]) -> dict[int, int]:
  counts: dict[int, int] = {}
  for pin in pin_map:
    key = int(pin)
    counts[key] = counts.get(key, 0) + 1
  return counts


def _problem_material_areas(
  pin_map: Sequence[int],
  *,
  widths: Sequence[float],
  pitch: float,
  radius: float,
  assigned: bool,
  moderator_material: int = 6,
) -> dict[int, float]:
  areas: dict[int, float] = {}
  for pin, count in _pin_counts(pin_map).items():
    if assigned:
      pin_areas = _pin_assigned_material_areas(
        pin,
        widths=widths,
        pitch=pitch,
        radius=radius,
        moderator_material=moderator_material,
      )
    else:
      pin_areas = _pin_true_material_areas(
        pin,
        pitch=pitch,
        radius=radius,
        moderator_material=moderator_material,
      )
    for material, area in pin_areas.items():
      areas[material] = areas.get(material, 0.0) + count * area
  return areas


def _material_area_corrections(
  target_areas: dict[int, float],
  assigned_areas: dict[int, float],
  *,
  total_area: float,
) -> tuple[MaterialAreaCorrection, ...]:
  entries = []
  for material in sorted(set(target_areas) | set(assigned_areas)):
    target = float(target_areas.get(material, 0.0))
    assigned = float(assigned_areas.get(material, 0.0))
    conservative = assigned > 0.0 or abs(target) <= 1.0e-14
    density_correction = target / assigned if assigned > 0.0 else None
    entries.append(
      MaterialAreaCorrection(
        material=material,
        material_name=PIN_NAMES.get(material, f"material {material}"),
        target_area=target,
        assigned_area=assigned,
        target_fraction=target / total_area if total_area > 0.0 else 0.0,
        assigned_fraction=assigned / total_area if total_area > 0.0 else 0.0,
        area_error=assigned - target,
        fraction_error=(assigned - target) / total_area if total_area > 0.0 else 0.0,
        density_correction=density_correction,
        conservative=conservative,
      )
    )
  return tuple(entries)


def center_sampled_mesh_recipe(
  pin_map: Sequence[int],
  *,
  pins_per_side: int,
  number: int = 7,
  pitch: float = PIN_PITCH,
  radius: float = PIN_RADIUS,
  correction_scope: str = "problem",
) -> CenterSampledMeshRecipe:
  """Build problem-level center-sampled mesh and density-correction metadata."""

  if len(pin_map) != pins_per_side * pins_per_side:
    raise ValueError("pin map size does not match pins_per_side")
  widths = uniform_pin_widths(number, pitch=pitch)
  edges = _edges_from_widths(widths)
  total_area = len(pin_map) * pitch * pitch
  target_areas = _problem_material_areas(
    pin_map,
    widths=widths,
    pitch=pitch,
    radius=radius,
    assigned=False,
  )
  assigned_areas = _problem_material_areas(
    pin_map,
    widths=widths,
    pitch=pitch,
    radius=radius,
    assigned=True,
  )
  return CenterSampledMeshRecipe(
    name=f"center_sampled_uniform_{int(number)}x{int(number)}_density_corrected",
    number=int(number),
    pitch=float(pitch),
    radius=float(radius),
    widths=tuple(widths),
    edges=tuple(edges),
    pins_per_side=int(pins_per_side),
    pin_count=len(pin_map),
    total_area=total_area,
    correction_scope=correction_scope,
    material_area_corrections=_material_area_corrections(
      target_areas,
      assigned_areas,
      total_area=total_area,
    ),
  )


def material_map_for_pin(
  central_material: int,
  *,
  moderator_material: int = 6,
  recipe: PinMeshRecipe | None = None,
) -> list[int]:
  """Return a row-major 7x7 material map for one C5G7 pin."""

  recipe = recipe or optimized_pin_recipe_7()
  if central_material == moderator_material:
    return [moderator_material] * (recipe.number * recipe.number)
  return [
    central_material if flag else moderator_material
    for flag in recipe.central_material_mask
  ]


def repeated_edges(widths: Sequence[float], count: int) -> list[float]:
  """Repeat one pin's widths over a row of pins."""

  edges = [0.0]
  for _ in range(int(count)):
    for width in widths:
      edges.append(edges[-1] + float(width))
  return edges


def expanded_materials(
  pin_map: Sequence[int],
  *,
  pins_per_side: int,
  recipe: PinMeshRecipe | None = None,
) -> list[int]:
  """Expand a pin map into a row-major fine-cell material map."""

  recipe = recipe or optimized_pin_recipe_7()
  if len(pin_map) != pins_per_side * pins_per_side:
    raise ValueError("pin map size does not match pins_per_side")

  pin_materials = {
    pin: material_map_for_pin(pin, recipe=recipe)
    for pin in set(int(value) for value in pin_map)
  }
  fine = []
  n = recipe.number
  for pin_j in range(pins_per_side):
    for sub_j in range(n):
      for pin_i in range(pins_per_side):
        pin = int(pin_map[pin_j * pins_per_side + pin_i])
        row = pin_materials[pin][sub_j * n:(sub_j + 1) * n]
        fine.extend(row)
  return fine


def center_sampled_expanded_materials(
  pin_map: Sequence[int],
  *,
  pins_per_side: int,
  recipe: CenterSampledMeshRecipe,
) -> list[int]:
  """Expand a pin map using center-sampled pin-cell material assignment."""

  if len(pin_map) != pins_per_side * pins_per_side:
    raise ValueError("pin map size does not match pins_per_side")

  pin_materials = {
    pin: center_sampled_material_map_for_pin(
      pin,
      widths=recipe.widths,
      pitch=recipe.pitch,
      radius=recipe.radius,
    )
    for pin in set(int(value) for value in pin_map)
  }
  fine = []
  n = recipe.number
  for pin_j in range(pins_per_side):
    for sub_j in range(n):
      for pin_i in range(pins_per_side):
        pin = int(pin_map[pin_j * pins_per_side + pin_i])
        row = pin_materials[pin][sub_j * n:(sub_j + 1) * n]
        fine.extend(row)
  return fine


def build_cartesian_mesh(dt, edges: Sequence[float], materials: Sequence[int]):
  """Build a Detran Mesh2D from fine edges and row-major materials."""

  return dt.CartesianMesh.from_coarse(
    x=list(edges),
    y=list(edges),
    materials=list(materials),
  )


def build_pin_mesh(dt, pin: int, recipe: PinMeshRecipe | None = None):
  recipe = recipe or optimized_pin_recipe_7()
  return build_cartesian_mesh(
    dt,
    recipe.edges,
    material_map_for_pin(pin, recipe=recipe),
  )


def build_lattice_mesh(
  dt,
  pin_map: Sequence[int],
  *,
  pins_per_side: int,
  recipe: PinMeshRecipe | None = None,
):
  recipe = recipe or optimized_pin_recipe_7()
  return build_cartesian_mesh(
    dt,
    repeated_edges(recipe.widths, pins_per_side),
    expanded_materials(pin_map, pins_per_side=pins_per_side, recipe=recipe),
  )


def build_center_sampled_problem_mesh(
  dt,
  pin_map: Sequence[int],
  *,
  pins_per_side: int,
  number: int = 7,
) -> tuple[object, CenterSampledMeshRecipe]:
  """Build a center-sampled Cartesian mesh and its correction recipe."""

  recipe = center_sampled_mesh_recipe(
    pin_map,
    pins_per_side=pins_per_side,
    number=number,
  )
  mesh = build_cartesian_mesh(
    dt,
    repeated_edges(recipe.widths, pins_per_side),
    center_sampled_expanded_materials(
      pin_map,
      pins_per_side=pins_per_side,
      recipe=recipe,
    ),
  )
  return mesh, recipe


def core_pin_map(
  assembly_map: Sequence[int],
  assembly_pin_maps: Sequence[Sequence[int]],
  *,
  assemblies_per_side: int = CORE_DIMENSION,
  pins_per_assembly_side: int = ASSEMBLY_DIMENSION,
) -> list[int]:
  """Expand an assembly map into the quarter-core pin map."""

  pins = []
  for assembly_j in range(assemblies_per_side):
    for pin_j in range(pins_per_assembly_side):
      for assembly_i in range(assemblies_per_side):
        assembly = int(assembly_map[assembly_j * assemblies_per_side + assembly_i])
        pin_map = assembly_pin_maps[assembly]
        offset = pin_j * pins_per_assembly_side
        pins.extend(pin_map[offset:offset + pins_per_assembly_side])
  return pins


def mesh_audit(recipe: PinMeshRecipe | None = None) -> dict:
  """Return the derivation quantities used by the report and tests."""

  recipe = recipe or optimized_pin_recipe_7()
  return recipe.to_dict()
