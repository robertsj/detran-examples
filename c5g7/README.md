# C5G7 Examples

This folder contains 2-D C5G7 benchmark inputs and small drivers for diffusion,
SN transport, native MOC, native CSG-MOC, CSG geometry audits, and MOC/SN
comparison runs.  The defaults favor short local checks over reference-quality
benchmark production runs.

Run examples from this directory with a current pydetran build on `PYTHONPATH`:

```bash
export PYTHONPATH=/home/robertsj/Research/detran/install/python:$PYTHONPATH
python3 run_pin_moc.py --index 0 --pin-mesh 7 --plot-style cells
python3 run_assembly_moc.py --index 0 --pin-mesh 7 --plot-style cells
python3 run_assembly_moc.py --equation ddmoc --index 0 --pin-mesh 7 --no-plots
python3 run_core_moc.py --pin-mesh 7 --plot-groups 0 --eigen-iters 300
python3 -m c5g7.run_case_csg_pin --track --solve-csg-moc --moc-spacing 2.0
python3 -m c5g7.run_case_sn_core --audit-only
```

The MOC drivers report `keff`, elapsed time, sweep count, and reference
eigenvalue differences when a matching reference value is available.  They can
write scalar-flux plots with `--plot-style cells` for compact matplotlib mesh
cell SVGs or `--plot-style tracks` for detailed track-segment SVGs.  Use
`--equation scmoc` for step characteristics or `--equation ddmoc` for
diamond-difference MOC.  Use `--plot-groups 0,2-4`, `--plot-groups all`,
`--plot-dir <path>`, `--plot-palette <name>`, or `--no-plots` to control plot
output.

Generated plot folders such as `moc_plots/` and `moc_plots_mesh7_cases/` are
local run artifacts and are ignored by git.  The timing CSVs and plot images in
those folders can be regenerated from the MOC drivers.

The CSG drivers build analytic circular-pin geometry with pin-level controls
for fuel radial divisions, moderator radial divisions, and azimuthal sectors.
By default, pins use one fuel radial zone, one moderator radial zone, and eight
azimuthal sectors with boundaries every 45 degrees.  The CSG drivers write JSON
and LaTeX row fragments under `csg_results/`; the self-contained report
`c5g7_csg_benchmark_report.tex` loads those fragments when present.
`--solve-csg-moc` runs native CSG MOC through the Detran eigenvalue manager.
With the default `--csg-overlay-mesh-variant none`, this uses one flat-source
region per analytic CSG region.  With `--csg-overlay-mesh-variant centered`,
`legacy`, or `optimized`, the true CSG tracks are split by Cartesian mesh planes
and each non-empty `(CSG region, overlay cell)` intersection becomes an FSR.
The overlay mesh contributes only FSR cuts; materials come from the analytic CSG
geometry.  `--solve-mesh-moc` remains available as a companion staircase mesh
comparison backend that solves directly on the SN-style mesh.

The SN/DD artifact drivers wrap the legacy mesh-backed transport path while
making the solver and mesh recipe explicit:

```bash
python3 -m c5g7.run_case_sn_pin --audit-only
python3 -m c5g7.run_case_sn_assembly --audit-only
python3 -m c5g7.run_case_sn_core --audit-only
```

Drop `--audit-only` to run transport.  The default SN backend is SLEPc GD
through the Arnoldi operator path with PETSc linear solves and MG-CMDSA
preconditioning where requested:

```bash
python3 -m c5g7.run_case_sn_core --eigen-backend slepc-gd --acceleration mgcmdsa
```

Fallback/comparison choices include `--eigen-backend gd`,
`--eigen-backend slepc-krylovschur`, `--eigen-backend cmfd`, and
`--eigen-backend power`, plus `--acceleration cmfd` or `--acceleration none`.
SN artifacts are written under `results/sn/`.

The default SN mesh variant is `--mesh-variant optimized`, currently a 7x7
automated volume-preserving recipe.  It uses the D4-symmetric 13-cell central
mask from the historical C5G7 staircase mesh, but computes the widths from a
one-dimensional overlap minimization subject to exact central-area preservation.
Use `--mesh-variant legacy` to call `PinCell.meshify(7, True)` instead.

## Source Files

| File | Purpose |
| --- | --- |
| `README.md` | This map of the C5G7 example folder. |
| `material_c5g7.py` | Builds the seven-group C5G7 material table with benchmark cross sections. |
| `pins_c5g7.py` | Builds the seven pin-cell definitions and fine meshes used by assemblies and pin cases. |
| `assemblies_c5g7.py` | Builds the UO2, MOX, and reflector assembly layouts from the pin definitions. |
| `core_c5g7.py` | Builds the quarter-core and larger repeated-core layouts from the assemblies. |
| `reference.py` | Stores C5G7 reference eigenvalues and pin-power regression data. |
| `csg_c5g7.py` | Builds true 2-D CSG pin, assembly, and quarter-core C5G7 geometry with radial/azimuthal pin subdivision. |
| `c5g7_csg_benchmark_report.tex` | Self-contained LaTeX benchmark report that inputs generated CSG result rows. |
| `plot_utils.py` | Shared mesh-cell and MOC track plotting helpers for scalar flux and mesh maps. |
| `moc_driver.py` | Shared native MOC problem construction, option parsing, solve execution, plotting, and summary reporting. |
| `run_moc.py` | Combined CLI for `--case pin`, `--case assembly`, or `--case core` native MOC runs. |
| `run_pin_moc.py` | Native MOC pin eigenproblem wrapper around `moc_driver.py`. |
| `run_assembly_moc.py` | Native MOC assembly eigenproblem wrapper around `moc_driver.py`. |
| `run_core_moc.py` | Native MOC quarter-core eigenproblem wrapper around `moc_driver.py`. |
| `run_case_csg_pin.py` | Generates CSG pin-cell geometry, tracking, and optional native CSG-MOC benchmark artifacts. |
| `run_case_csg_assembly.py` | Generates CSG assembly geometry, tracking, and optional native CSG-MOC benchmark artifacts. |
| `run_case_csg_core.py` | Generates CSG quarter-core geometry, tracking, and optional native CSG-MOC benchmark artifacts. |
| `sn_mesh.py` | Generates and audits the automated 7x7 volume-preserving SN pin mesh. |
| `sn_driver.py` | Shared native SN/DD artifact runner with PETSc/SLEPc and acceleration options. |
| `run_case_sn_pin.py` | Generates SN/DD pin-cell benchmark artifacts. |
| `run_case_sn_assembly.py` | Generates SN/DD assembly benchmark artifacts. |
| `run_case_sn_core.py` | Generates SN/DD quarter-core benchmark artifacts. |
| `run_core_compare.py` | Solves the quarter-core with native MOC and SN on the same mesh and reports normalized flux differences. |
| `run_pin.py` | Original pin eigenproblem driver using the transport/eigenvalue input helpers. |
| `run_pin_diffusion.py` | Diffusion pin eigenproblem smoke driver. |
| `run_assembly.py` | Original reflected assembly eigenproblem driver using DD SN transport. |
| `run_assembly_diffusion.py` | Diffusion assembly eigenproblem smoke driver with group-0 plotting. |
| `run_assembly_fixed.py` | Fixed-source assembly experiment for multigroup preconditioner behavior. |
| `run_core.py` | Original DD SN quarter-core eigenproblem with MG-CMDSA-style solver settings. |
| `run_core_diffusion.py` | Diffusion quarter-core eigenproblem smoke driver with group-0 plotting. |
| `run_large_core.py` | Larger repeated-core DD SN eigenproblem stress case. |

## Notes

Pin and assembly MOC cases use reflected boundaries on all sides.  The core MOC
case uses the same quarter-core convention as `run_core.py`: west and south
reflected, east and north vacuum.

The mesh density is controlled by `--pin-mesh`; use `--pin-mesh 7` for the
current C5G7 mesh-7 variant.  SN runs also accept `--mesh-variant centered`,
which assigns each pin subcell by its center point and emits a problem-level
material-area correction table.  MOC runs accept the same centered variant so
MOC and SN can be compared on the same refinement ladder.  The full core can
need a higher power-iteration cap than the short default, so `--eigen-iters
300` is a better convergence timing run setting for the mesh-7 core.
