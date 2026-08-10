# C5G7 Examples

This folder contains 2-D C5G7 benchmark inputs and small drivers for diffusion,
SN transport, native MOC, and MOC/SN comparison runs.  The
defaults favor short local checks over reference-quality benchmark production
runs.

Run examples from this directory with a current pydetran build on `PYTHONPATH`:

```bash
export PYTHONPATH=/home/robertsj/Research/detran/install/python:$PYTHONPATH
python3 run_pin_moc.py --index 0 --pin-mesh 7 --plot-style cells
python3 run_assembly_moc.py --index 0 --pin-mesh 7 --plot-style cells
python3 run_assembly_moc.py --equation ddmoc --index 0 --pin-mesh 7 --no-plots
python3 run_core_moc.py --pin-mesh 7 --plot-groups 0 --eigen-iters 300
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

## Source Files

| File | Purpose |
| --- | --- |
| `README.md` | This map of the C5G7 example folder. |
| `material_c5g7.py` | Builds the seven-group C5G7 material table with benchmark cross sections. |
| `pins_c5g7.py` | Builds the seven pin-cell definitions and fine meshes used by assemblies and pin cases. |
| `assemblies_c5g7.py` | Builds the UO2, MOX, and reflector assembly layouts from the pin definitions. |
| `core_c5g7.py` | Builds the quarter-core and larger repeated-core layouts from the assemblies. |
| `reference.py` | Stores C5G7 reference eigenvalues and pin-power regression data. |
| `plot_utils.py` | Shared mesh-cell and MOC track plotting helpers for scalar flux and mesh maps. |
| `moc_driver.py` | Shared native MOC problem construction, option parsing, solve execution, plotting, and summary reporting. |
| `run_moc.py` | Combined CLI for `--case pin`, `--case assembly`, or `--case core` native MOC runs. |
| `run_pin_moc.py` | Native MOC pin eigenproblem wrapper around `moc_driver.py`. |
| `run_assembly_moc.py` | Native MOC assembly eigenproblem wrapper around `moc_driver.py`. |
| `run_core_moc.py` | Native MOC quarter-core eigenproblem wrapper around `moc_driver.py`. |
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

The mesh variant is controlled by `--pin-mesh`; use `--pin-mesh 7` for the
current C5G7 mesh-7 variant.  The full core can need a higher power-iteration
cap than the short default, so `--eigen-iters 300` is a better convergence
timing run setting for the mesh-7 core.
