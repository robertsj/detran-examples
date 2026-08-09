# Multigroup Diffusion Preconditioners for Multiplying Fixed-Source Problems

These examples reproduce, to the degree practical with the current built-in
solvers, results from:

    J. Roberts and B. Forget, "Multigroup Diffusion Preconditioners for
    Multiplying, Fixed-Source Transport Problems," Journal of Computational
    Physics 274, 1 October 2014, Pages 455-472.

Paper DOI: <https://doi.org/10.1016/j.jcp.2014.06.034>

The original paper used PETSc for some linear solver and preconditioner pieces.
The default archival workflow deliberately uses Detran's built-in `callow`
solvers only: GMRES plus built-in ILU0/ILUT/Jacobi-style preconditioning where
applicable.  A PETSc comparison path is also available when detran is built with
PETSc support.

The C5G7 setup uses the paper's uniform fast-group fixed boundary source on the
west face.  The built-in backend uses a full-fill ILUT solve for the paper rows
labeled LU, because callow does not currently provide a direct LU
preconditioner; PETSc mode maps those rows to PETSc LU.

Run a quick smoke solve:

```bash
python3 run_all.py
```

Run the paper reproduction smoke.  This runs a cheap C5G7 fixed-boundary
transport solve, assembles reduced BWR operators for all preconditioner kinds,
and checks that the full BWR paper geometry still gives the published
3267-by-3267 operator size:

```bash
python3 run_all.py reproduce-paper --paper-mode smoke
```

Run the full archival workflow:

```bash
python3 run_all.py all --quadrature 8 --tolerance 1e-6
```

Run a current reproduction for one paper tolerance:

```bash
python3 run_all.py reproduce-paper --paper-mode current --quadrature 8 --tolerance 1e-6
```

Run the full C5G7 paper table workflow for all three published tolerances:

```bash
python3 run_all.py reproduce-paper --paper-mode full --quadrature 8
```

Run the PETSc-backed comparison workflow:

```bash
python3 run_all.py all --linear-backend petsc --quadrature 8 --tolerance 1e-6
```

To compare PETSc timings against the already-calibrated built-in tolerances,
skip PETSc recalibration and run only the case timings:

```bash
python3 run_all.py run --linear-backend petsc --calibration-source builtin --quadrature 8 --tolerance 1e-6
python3 run_all.py table --linear-backend petsc --calibration-source builtin --quadrature 8 --tolerance 1e-6
```

The full workflow writes local pickle/text artifacts under this directory:
reference flux and pin-power data, calibrated tolerances, timing/sweep data, and
a LaTeX table summarizing the selected solver/preconditioner cases.  PETSc
artifacts use a `_petsc` suffix so the built-in archival data stays intact.
Reference and calibration files include the quadrature and tolerance in their
names so all three paper tolerances can coexist.

## Updated local results

The `updated_results/` subdirectory contains a compact archive of the local
run data used for the current example report.  It intentionally stores only the
small plotting inputs: C5G7 timing/sweep rows, BWR spectral diagnostics, and
BWR eigenvalue caches.  Dense BWR operator matrices and raw logs are generated
artifacts and are ignored by git.

Rebuild the PDF report from the archived plotting data with:

```bash
python3 run_all.py updated-results
```

Refresh the compact archive from local generated run artifacts, then rebuild
the report, with:

```bash
python3 run_all.py archive-results
python3 run_all.py updated-results
```

The committed PDF is `updated_results/updated_results.pdf`.

## Published paper data

The recovered manuscript tables and the small BWR spectral benchmark data are
stored in `paper_results.json`.  Print the published C5G7 table for a selected
tolerance with:

```bash
python3 run_all.py paper-table --tolerance 1e-6
```

After generating current timing/sweep data with `run_all.py run`, compare the
current built-in results to the paper table:

```bash
python3 run_all.py run --quadrature 8 --tolerance 1e-6
python3 run_all.py paper-compare --quadrature 8 --tolerance 1e-6
```

Long C5G7 runs checkpoint after each completed row.  If a row is interrupted,
resume without losing completed rows:

```bash
python3 run_all.py run --quadrature 8 --tolerance 1e-4 --resume
```

To isolate known slow rows or fill only missing rows, select explicit case keys:

```bash
python3 run_all.py run --quadrature 8 --tolerance 1e-4 --cases dsa_lu dsa_ilu --resume
python3 run_all.py paper-compare --quadrature 8 --tolerance 1e-4 --cases gmres cmdsa cmdsa_s
```

Calibration artifacts are tolerance-specific, so the three paper calibrations
can be launched in separate shells before doing timing-sensitive `run` commands:

```bash
python3 run_all.py calibrate --quadrature 8 --tolerance 1e-4
python3 run_all.py calibrate --quadrature 8 --tolerance 1e-6
python3 run_all.py calibrate --quadrature 8 --tolerance 1e-8
```

For recovered paper tolerances, `--calibration-seed paper` starts from the
published residuals and verifies them with actual transport solves before
accepting the calibration checkpoint:

```bash
python3 run_all.py calibrate --linear-backend petsc --calibration-seed paper --quadrature 8 --tolerance 1e-8 --resume
```

Print the published BWR spectral table with:

```bash
python3 run_all.py paper-bwr
```

## BWR transport operators

The BWR operator path reconstructs the three-group benchmark from the recovered
paper data.  The default is a small, finite diagnostic case: the central pin
plus moderator ring, with the same 3-by-3 fine subdivision per pin used by the
paper's full-size setup.

```bash
python3 run_all.py bwr-operator --bwr-eigs --bwr-gmres
python3 run_all.py bwr-operator --bwr-operator-kind mgdsa --bwr-eigs --bwr-gmres
python3 run_all.py bwr-operator --bwr-operator-kind cmdsa --bwr-eigs --bwr-gmres
python3 run_all.py bwr-operator --bwr-operator-kind cmdsa-s --bwr-eigs --bwr-gmres
python3 run_all.py bwr-operator --bwr-operator-kind tcmdsa --bwr-eigs --bwr-gmres
```

Use `--bwr-paper-size` to select the full paper geometry and quadrature.  In
size-only mode this should report the published 3267-by-3267 operator:

```bash
python3 run_all.py bwr-operator --bwr-paper-size --bwr-size-only
```

Full paper-size dense assembly is intentionally opt-in.  It writes a compressed
`.npz` containing the dense operator, JSON summary, and coarse material map:

```bash
python3 run_all.py bwr-operator --bwr-paper-size --bwr-operator-kind transport
```

The orchestrated full-paper command only checks the full BWR size unless this
flag is added:

```bash
python3 run_all.py reproduce-paper --paper-mode full --include-full-bwr-dense
```

The preconditioned operator kinds are `mgdsa`, `cmdsa`, `cmdsa-s`, and
`tcmdsa`.  The built-in diffusion solve uses no additional preconditioner by
default; `--bwr-diffusion-pc ilu0` or `--bwr-diffusion-pc ilut` can be used for
experiments.  For diagnostics, `--bwr-disable-pc-fission` omits fission from
the CMDSA-family BWR preconditioner while leaving the transport operator setup
otherwise unchanged.
