"""Plotting helpers shared by the C5G7 example drivers."""

from pathlib import Path

import numpy as np
try:
  import matplotlib
  import matplotlib.pyplot as plt
except ImportError:
  matplotlib = None
  plt = None

try:
  from detran_viz import track_region_quantity_to_svg
except ImportError:
  try:
    from detran import track_region_quantity_to_svg
  except ImportError:
    track_region_quantity_to_svg = None

try:
  from detran import mesh_track_product
except ImportError:
  mesh_track_product = None

def has_matplotlib():
  return plt is not None and matplotlib is not None

def require_matplotlib():
  if plt is None or matplotlib is None:
    raise RuntimeError("matplotlib is required for C5G7 plotting helpers")

def has_track_quantity_svg():
  return (
    mesh_track_product is not None
    and track_region_quantity_to_svg is not None
  )

def require_track_quantity_svg():
  if not has_track_quantity_svg():
    raise RuntimeError(
      "C5G7 track-segment flux plots require detran.mesh_track_product and "
      "detran_viz.track_region_quantity_to_svg; rebuild/install current pydetran"
    )

def mesh_quantity_grid(obj, values):
  """Return a fine-mesh scalar grid in Detran's native cell ordering."""
  values = np.asarray(values, dtype=float)
  if values.size != obj.number_cells():
    raise ValueError(
      "mesh quantity has {0} cells, expected {1}".format(
        values.size,
        obj.number_cells(),
      )
  )
  return values.reshape(obj.number_cells_y(), obj.number_cells_x())

def mesh_cell_index_grid(obj):
  """Return a fine-mesh cell id grid in Detran's native cell ordering."""
  return np.arange(obj.number_cells(), dtype=np.int64).reshape(
    obj.number_cells_y(),
    obj.number_cells_x(),
  )

def write_moc_flux_svgs(
  obj,
  result,
  output_dir,
  prefix,
  groups=None,
  width=900,
  palette="viridis",
  style="tracks",
):
  """Write one scalar-flux SVG per requested group for a MOC result."""
  style = str(style).lower()
  if style not in ("tracks", "cells"):
    raise ValueError("style must be 'tracks' or 'cells'")
  if style == "tracks":
    require_track_quantity_svg()
  else:
    require_matplotlib()
  output_dir = Path(output_dir)
  output_dir.mkdir(parents=True, exist_ok=True)

  if groups is None:
    group_count = getattr(getattr(result, "problem", None), "groups", None)
    if group_count is None:
      raise ValueError("groups must be supplied when result.problem.groups is absent")
    groups = range(int(group_count))

  if style == "tracks":
    product = mesh_track_product(obj, result.solver.input(), result.solver.quadrature())
    index = "region_id"
    filename_template = "{0}_track_flux_g{1}.svg"
  else:
    product = None
    index = None
    filename_template = "{0}_cell_flux_g{1}.svg"

  plot_files = []
  for group in groups:
    group = int(group)
    flux = np.asarray(result.flux(group), dtype=float)
    if flux.shape[0] != obj.number_cells():
      raise ValueError(
        "group {0} flux has {1} cells, expected {2}".format(
          group,
          flux.shape[0],
          obj.number_cells(),
        )
    )
    filename = filename_template.format(prefix, group)
    path = output_dir / filename
    title = "{0} scalar flux group {1}".format(prefix, group)
    label = "scalar flux g{0}".format(group)
    if style == "tracks":
      svg = track_region_quantity_to_svg(
        product,
        flux,
        index=index,
        width=width,
        palette=palette,
        colorbar_label=label,
        title=title,
      )
      path.write_text(svg)
    else:
      save_mesh_quantity_svg(
        obj,
        flux,
        path,
        width=width,
        cmap=palette,
        title=title,
        colorbar_label=label,
      )
    plot_files.append(str(path))
  return plot_files

def plot_grid_quantity(
  x_edges,
  y_edges,
  quantity,
  *,
  ax=None,
  cmap="viridis",
  title=None,
  colorbar=True,
  colorbar_label=None,
  edges=False,
):
  """Plot a cell-centered scalar grid using explicit mesh edges."""
  require_matplotlib()
  x_edges = np.asarray(x_edges, dtype=float)
  y_edges = np.asarray(y_edges, dtype=float)
  quantity = np.asarray(quantity, dtype=float)
  expected = (y_edges.size - 1, x_edges.size - 1)
  if quantity.shape != expected:
    raise ValueError(
      "quantity grid shape is {0}, expected {1}".format(quantity.shape, expected)
    )

  if ax is None:
    fig, ax = plt.subplots()
  else:
    fig = ax.figure

  if edges:
    edgecolors = "black"
    linewidth = 0.15
  else:
    edgecolors = "none"
    linewidth = 0.0

  image = ax.pcolormesh(
    x_edges,
    y_edges,
    np.ma.masked_invalid(quantity),
    cmap=cmap,
    shading="flat",
    edgecolors=edgecolors,
    linewidth=linewidth,
  )
  ax.set_aspect("equal")
  ax.set_xlabel("x")
  ax.set_ylabel("y")
  if title:
    ax.set_title(title)
  if colorbar:
    bar = fig.colorbar(image, ax=ax)
    if colorbar_label:
      bar.set_label(colorbar_label)
  return fig, ax, image

def save_grid_quantity_svg(
  x_edges,
  y_edges,
  quantity,
  path,
  *,
  width=900,
  dpi=100,
  cmap="viridis",
  title=None,
  colorbar_label=None,
  edges=False,
):
  """Save a cell-centered scalar grid as a matplotlib SVG with colorbar."""
  require_matplotlib()
  x_edges = np.asarray(x_edges, dtype=float)
  y_edges = np.asarray(y_edges, dtype=float)
  span_x = max(float(x_edges[-1] - x_edges[0]), 1.0e-12)
  span_y = max(float(y_edges[-1] - y_edges[0]), 1.0e-12)
  width_in = max(float(width), 100.0) / float(dpi)
  height_in = max(2.5, width_in * span_y / span_x)

  with matplotlib.rc_context({"svg.fonttype": "none"}):
    fig, ax = plt.subplots(figsize=(width_in, height_in), dpi=dpi)
    try:
      plot_grid_quantity(
        x_edges,
        y_edges,
        quantity,
        ax=ax,
        cmap=cmap,
        title=title,
        colorbar_label=colorbar_label,
        edges=edges,
      )
      fig.tight_layout()
      fig.savefig(str(path), format="svg", bbox_inches="tight")
    finally:
      plt.close(fig)
  return str(path)

def save_mesh_quantity_svg(obj, values, path, **kwargs):
  """Save a Detran fine-mesh scalar quantity as a matplotlib SVG."""
  x, y = mesh_axes(obj)
  return save_grid_quantity_svg(x, y, mesh_quantity_grid(obj, values), path, **kwargs)

def plot_mesh_map(obj, map_key, edges=False):
  """Plot an integer mesh map."""
  if not has_matplotlib():
    print("matplotlib is not available; skipping plot.")
    return
  m_map = np.asarray(obj.mesh_map(map_key), dtype=float)
  if m_map.size == 0:
    print("Oops! Map has size zero!")
    return
  x, y = mesh_axes(obj)
  fig, ax, image = plot_grid_quantity(
    x,
    y,
    m_map.reshape(obj.number_cells_y(), obj.number_cells_x()),
    cmap="tab20",
    title=map_key,
    colorbar_label=map_key,
    edges=edges,
  )
  plt.show()
  return fig, ax, image


def plot_flux(obj, f):
  """Plot a flux or flux-shaped vector on the mesh."""
  if not has_matplotlib():
    print("matplotlib is not available; skipping plot.")
    return
  x, y = mesh_axes(obj)
  fig, ax, image = plot_grid_quantity(
    x,
    y,
    mesh_quantity_grid(obj, f),
    cmap="hot",
    colorbar_label="flux",
  )
  plt.show()
  return fig, ax, image

#  plt.imshow(f, interpolation='nearest', cmap=plt.cm.hot)
#  plt.title('flux map')
#  plt.colorbar()
#  plt.show()

def mesh_axes(obj):
  """Return the fine-mesh edges defining cells for plotting."""
  x = np.zeros(obj.number_cells_x()+1)
  y = np.zeros(obj.number_cells_y()+1)
  for i in range(0, obj.number_cells_x()):
    x[i + 1] = x[i] + obj.dx(i)
  for j in range(0, obj.number_cells_y()):
    y[j + 1] = y[j] + obj.dy(j)
  return (x, y)
