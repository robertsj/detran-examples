# CSG Examples

This folder starts a CSG inspection hierarchy for Detran.

The first model, `pwr17_assembly.py`, is a VERA-style 17x17 PWR assembly
scaffold.  It uses public VERA Core Physics Benchmark progression dimensions
for the fuel rod, lattice, spacer grids, assembly height, and RCCA control rod
data, but it is not yet a certified transcription of a benchmark input deck.

The `spheroid_performance.py` playground uses the live pydetran CSG and
tracking API rather than the older sampled-ray scaffold.  It embeds analytic
spheres and oblong spheroids in the unit box, scales their semi-axes to stay in
the `[0.01, 0.99]` guard band, and reports separate build, slice, and tracking
timings along with track-derived volume estimates.  When an output directory is
requested, it also writes exact analytic conic SVGs for the midplane slice and
orthographic projection of the spheroid, so sampled products can be compared to
native C++ quadratic feature boundaries and hand-calculable geometry.

The hierarchy has deliberate exit ramps:

- `homogeneous_model()` gives a one-material assembly-sized box.
- `fuel_pin_2d_model()` gives a single 2-D fuel pin cell.
- `assembly_2d_model(control_position=...)` gives a 17x17 x-y slice.
- `assembly_3d_model(control_position=..., include_spacer_grids=...)` adds
  axial height, spacer-grid bands, nozzles, and variable control insertion.

The control convention is `control_position=0.0` for withdrawn and `1.0` for
fully inserted.  Intermediate values move the poison stack using the VERA RCCA
step size and maximum step count.

Run the example from this directory with pydetran on `PYTHONPATH`:

```bash
PYTHONPATH=/home/robertsj/Research/detran/pydetran \
python3 pwr17_assembly.py --output-dir out --control-position 1.0
```

Run the spheroid performance playground against a locally built pydetran:

```bash
PYTHONPATH=/home/robertsj/Research/detran/install/python \
LD_LIBRARY_PATH=/home/robertsj/Research/detran/install/lib \
python3 spheroid_performance.py --case all --spacing 0.05 --output-dir out
```

The spheroid playground writes four SVG inspection products per case when
``--output-dir`` is set:

- ``*_slice_xy.svg``: sampled region map from the native CSG geometry.
- ``*_tracks_xy.svg``: deterministic C++ tracker segments.
- ``*_analytic_slice_xy.svg``: exact conic for the requested slice plane.
- ``*_analytic_projection_xy.svg``: exact orthographic projection footprint.

The analytic products are geometry features, not a reduced transport model.
They are useful for visual checks, coordinate indexing tests, and mesh/AMR
seeding. A true flattened 2-D tracking model would need an explicit material
collapse policy before it should be used as a physics shortcut.

The PWR17 scaffold outputs deterministic inspection products:

- PNG slices colored by semantic region/material.
- SVG ray traces colored by the region crossed by each ray segment, with one
  file per sampled angle and axial plane.
- Optional HDF5 at `out/pwr17_csg_products.h5` with `/geometry/source`,
  `/geometry/regions`, `/geometry/products/slices`, and
  `/geometry/products/rays`.

By default, generated assembly pixel maps use `1445 x 1445` cells, five times
the original scaffold resolution of `289 x 289`.  The smaller homogeneous and
pin-cell exit ramps are also scaled by five.  Use `--nx` and `--ny` for cheaper
smoke runs.

The default ray angles are `0`, `12`, `30`, `45`, `60`, and `78` degrees.  Use
`--ray-angle` or `--ray-angles` followed by one or more values to override the
set.

Fixed-grid slices are PNG by default because a high-resolution pixel map is a
raster product.  Use `--slice-format svg` or `--slice-format both` for debug
SVGs where each rectangle keeps region metadata and can be inspected or
overlaid without blur.  The ray SVGs are different: they are genuinely vector
products, with one line element per sampled region segment.

The PWR17 scaffold intentionally emits gridded products rather than a full
Detran C++ `Geometry`.  The spheroid playground is the smaller live-backend
counterpart: it uses pydetran CSG geometry and the C++ trackers directly while
keeping the same reusable product vocabulary of slices and tracks.

## Notebook walkthroughs

Two Jupyter notebooks provide ready-to-run visual walkthroughs:

- `spheroid_viz_walkthrough.ipynb` builds the live pydetran spheroid CSG
  cases, compares sampled slices against exact analytic conics, displays
  projected tracker segments, and runs small performance/volume checks.
- `pwr17_assembly_viz_walkthrough.ipynb` walks through the PWR scaffold from a
  pin-cell slice to 2-D assembly slices, 3-D axial planes, sampled ray traces,
  and a compact SVG output bundle.

Run them with the installed pydetran runtime visible to the kernel:

```bash
export DETRAN_ROOT=/home/robertsj/Research/detran
export PATH=$DETRAN_ROOT/install/bin:$PATH
export LD_LIBRARY_PATH=$DETRAN_ROOT/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$DETRAN_ROOT/install/python:$PYTHONPATH
cd $DETRAN_ROOT/detran-examples/csg
jupyter notebook
```

Notebook-generated files are written below `notebook_outputs/`, keeping them
separate from the command-line `out/` examples.
