# C5G7 2-D Benchmark Report Plan

This plan turns the current C5G7 example work into a reproducible benchmark
report program.  The immediate scope is 2-D.  The 3-D estimates and tracking
work remain relevant background, but they should not block the 2-D report.

## Objectives

- Produce a self-contained LaTeX report for the C5G7 benchmark family.
- Compare Detran native CSG transport, native MOC on mesh geometry, and native
  SN/DD transport on the existing 7x7 volume-preserving pin mesh.
- Preserve enough run metadata that every table entry can be regenerated from
  small `run_case_*.py` drivers.
- Separate short smoke/default cases from reference-quality and refinement
  studies.
- Prepare the report to include Ganapol-style sequence acceleration once
  refinement sequences are available.

## Report Outline

1. Benchmark definition
   - Summarize the C5G7 problem variants used here: pin cells, reflected
     assemblies, and the 2-D quarter-core.
   - Document materials, pin radius, pin pitch, assembly pitch, core layout,
     and boundary conditions.
   - State reference quantities: eigenvalues, pin powers, normalization rules,
     and uncertainty/source status for published Monte Carlo reference data.
   - Distinguish official C5G7 reference data from Detran regression values and
     other external deterministic results.

2. Detran models and discretizations
   - Native CSG geometry: analytic circular-pin CSG with controls
     `fuel_radial_divisions`, `moderator_radial_divisions`, and
     `azimuthal_slices`.
   - Native CSG solve: direct region-track MOC through the CSG solver-manager
     bridge.  The first bridge uses one synthetic mesh cell per CSG region and
     a CSG-generated `TrackCollection`; exact analytic region volumes and
     direct pin-power reductions remain refinement work.
   - Native MOC mesh solve: current mesh-backed native MOC drivers using
     `pin_mesh`, track spacing, polar/azimuth quadrature, and equation choice.
   - Native SN/DD solve: the existing 7x7 volume-preserving pin mesh used by
     `run_pin.py`, `run_assembly.py`, and `run_core.py`, with the documented
     DD/SN angular quadratures and solver settings.
   - Mesh generation: define both the legacy `PinCell.meshify(number, flag)`
     path and the proposed optimized volume-preserving mesh path.

3. Results
   - Tables for representative and best-bet calculations:
     `k_eff`, pcm error, pin-power RMS error, pin-power max error, iterations,
     sweeps, wall time, thread count, and memory when available.
   - Separate tables for pin, assembly, and core cases.
   - Include Detran native CSG, Detran native MOC, Detran native SN/DD, and
     selected external published results where citation-quality data are
     available.
   - Add refinement tables over geometry, mesh, track spacing, angular
     quadrature, and solver tolerances.
   - Add sequence-accelerated estimates only when the raw sequence is clean
     enough to justify extrapolation.

4. Conclusions
   - Summarize best current Detran estimates and errors.
   - Identify dominant error sources: geometry, angular quadrature, track
     spacing, mesh staircase error, and solver tolerances.
   - State which model is recommended for production verification and which
     model is recommended for reference-quality studies.

## Run Artifact Contract

Every report row should come from a JSON artifact and a LaTeX row fragment under
`c5g7/results/` or a method-specific subdirectory.

Each JSON result should include at least:

- `case`: `pin`, `assembly`, or `core`
- `index`: pin or assembly index when applicable
- `method`: `native_csg`, `native_moc`, `native_sn`, or `external`
- `geometry_model`: concise geometry identifier
- `discretization`: structured mesh/CSG/track/quadrature settings
- `solver`: equation, iteration method, tolerances, limits, acceleration
- `reference`: reference value name, value, source, and uncertainty if known
- `results`: eigenvalue, pin powers, errors, iteration counts, timing, memory
- `artifacts`: paths to row fragments, plots, logs, and raw stdout when kept

## Implementation Work Packages

### 1. Stabilize build and result directories

- Create one active pydetran build configured against
  `/home/robertsj/Research/detran/pydetran`, not the older worktree.
- Keep generated results ignored by git unless they are intentionally curated.
- Use one canonical result root, probably `c5g7/results/`, with subdirectories
  `csg/`, `moc/`, `sn/`, `external/`, and `figures/`.

### 2. Simplify geometry helpers

The helper layer should stay declarative and small:

- `PinCellSpec`: pitch, central radius, central material, outer material, name.
- `AssemblySpec`: dimension, pin specifications, pin map, name.
- `CoreSpec`: dimension, assembly specifications, assembly map, boundary.

Each spec can emit:

- Native CSG geometry.
- Mesh geometry using a selected pin mesh recipe.
- Metadata tables for report auditing.

This can be implemented either by revising existing `PinCell`, `Assembly`, and
`Core` helpers or by adding a thin Python-only C5G7 scaffold.  Prefer the
smallest change that gives shared definitions for CSG, MOC, and SN.

### 3. Add native SN/DD run-case drivers

Create modern artifact-writing wrappers around the existing legacy scripts:

- `run_case_sn_pin.py`
- `run_case_sn_assembly.py`
- `run_case_sn_core.py`

Initial settings should reproduce the documented legacy runs:

- Pin: 7x7 volume-preserving mesh, DD/SN or diffusion as selected, QR-style
  quadrature from the existing driver comments.
- Assembly: 7x7 volume-preserving mesh, DD, reflected boundaries, quadrature of
  order 50 as currently documented.
- Core: 7x7 volume-preserving mesh, DD, west/south reflected and east/north
  vacuum, ASQR-ASDR polar 3 / azimuth 6, MGCMDSA settings from `run_core.py`.

The wrappers should write JSON, LaTeX rows, and any pin-power/error maps.

### 4. Harden native CSG solve bridge

The current CSG code builds analytic geometry, tracks it, and can perform
multigroup MOC eigenvalue solves by feeding CSG region tracks into the Detran
manager stack.  The first implementation is intentionally conservative: it
constructs a region-indexed synthetic mesh with one cell per CSG region and
normalizes CSG track-product region volumes to the analytic domain area.

Remaining hardening pieces:

- Exact analytic CSG region volumes where possible, or a documented high-order
  deterministic volume estimator.
- Region or pin-power tally reduction from CSG segments.
- Refinement studies on C5G7 pin, assembly, and core cases.
- Boundary and regression checks for reflected pin/assembly cases and the
  quarter-core reflect/vacuum boundary mix.

The report rows can now contain native CSG-MOC `k_eff` values for cases run
with `--solve-csg-moc`; rows without that flag remain geometry/tracking audits.

### 5. Investigate optimized volume-preserving NxN pin meshes

The legacy 7x7 mesh is useful and should be reported.  A better automated mesh
generator may be worth adding as an optional path.

Near-term baseline now available:

- `--mesh-variant centered` assigns each Cartesian pin subcell by its center
  point.
- The driver compares assigned material areas against known analytic C5G7
  material areas for the selected pin, assembly, or core.
- For each material with nonzero assigned area, the generated problem-specific
  material table scales macroscopic cross sections by
  `true_area / assigned_area`.
- If a material has true area but no assigned cells, the audit marks the recipe
  non-conservative because no density correction can recover a missing material.

Candidate algorithm:

1. Restrict to symmetric Cartesian grids for a square pin cell.
2. For fixed `N`, parameterize positive interval widths rather than raw grid
   lines, enforcing total pitch exactly.
3. Enumerate or choose a symmetric material mask for fuel/central material cells.
4. Compute exact rectangle-circle overlap areas, or use a high-order deterministic
   quadrature fallback for each candidate cell.
5. Minimize an objective:
   - exact or near-exact preservation of fuel and moderator area,
   - optional preservation of annular/radial CSG areas,
   - maximum cell-size uniformity,
   - symmetry and monotonicity of grid widths,
   - small material-boundary error against the true circle.
6. Emit the selected grid and material map as a reusable mesh recipe.

This is essentially a constrained optimization problem:

`minimize cell_size_variance + boundary_error + penalties`

subject to:

- total pitch preserved,
- target material areas preserved,
- fixed `N x N` cells,
- valid positive grid intervals.

Start with `N = 7` to compare directly against the legacy mesh, then generate
`N = 3, 5, 7, 9, 11` refinement ladders if the method behaves smoothly.

### 6. External references and protocols

- Locate and cite the official OECD/NEA C5G7 benchmark report.
- Record MCNP reference eigenvalues and pin-power data with uncertainties where
  available.
- Add a small `external_results.json` or BibTeX-backed table for published
  deterministic MOC/SN comparisons.
- Follow common presentation style:
  - eigenvalue errors in pcm,
  - normalized pin powers,
  - RMS and maximum relative pin-power errors,
  - clear geometry/angle/tolerance metadata in every table.

### 7. Refinement and sequence acceleration

Initial refinement axes:

- CSG radial divisions in fuel and moderator.
- CSG azimuthal sectors.
- MOC track spacing.
- MOC polar and azimuthal quadrature.
- SN mesh size and angular quadrature.
- Solver tolerances.

Sequence acceleration should be added after raw sequences exist:

- Store every raw sequence.
- Apply Aitken, Richardson/Romberg-style, Wynn epsilon, or Ganapol-style
  acceleration only when monotone/asymptotic behavior is plausible.
- Report both the raw best calculation and the accelerated estimate.
- Flag extrapolated values clearly as extrapolated, not direct solves.

## Suggested Execution Order

1. Rebuild pydetran against the active edited checkout and rerun focused tests.
2. Move the current CSG artifact paths to the canonical results layout.
3. Add SN/DD artifact-writing wrappers around the existing 7x7 drivers.
4. Expand the LaTeX report skeleton to the full four-section benchmark report.
5. Add the optimized 7x7 mesh prototype and compare it to legacy `meshify(7, True)`.
6. Harden the CSG solve bridge on C5G7 pin, assembly, and core refinement runs.
7. Add refinement driver generation and sequence-acceleration postprocessing.
8. Curate external references and finish the benchmark narrative.

## Current Known State

- The CSG scaffold currently produces geometry/tracking audit artifacts for pin,
  assembly, and core cases, and can produce native CSG-MOC transport artifacts
  when run with `--solve-csg-moc`.
- Native CSG-MOC eigenvalue transport has smoke coverage through pydetran and
  the C5G7 pin driver; reference-quality C5G7 production rows still need
  refinement settings, exact/validated volume treatment, and pin-power tallies.
- The legacy 7x7 SN/DD path exists in older drivers and should be wrapped rather
  than rewritten first.
- The current LaTeX report compiles, but it is still a short scaffold and needs
  the full benchmark sections above.
