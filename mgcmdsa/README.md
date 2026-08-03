# Multigroup Diffusion Preconditioners for Multiplying Fixed-Source Problems

These examples reproduce, to the degree practical with the current built-in
solvers, results from:

    J. Roberts and B. Forget, "Multigroup Diffusion Preconditioners for
    Multiplying, Fixed-Source Transport Problems," Journal of Computational
    Physics 274, 1 October 2014, Pages 455-472.

The original paper used PETSc for some linear solver and preconditioner pieces.
The default archival workflow deliberately uses Detran's built-in `callow`
solvers only: GMRES plus built-in ILU0/ILUT/Jacobi-style preconditioning where
applicable.  A PETSc comparison path is also available when detran is built with
PETSc support.

Run a quick smoke solve:

```bash
python3 run_all.py
```

Run the full archival workflow:

```bash
python3 run_all.py all --quadrature 8 --tolerance 1e-6
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
