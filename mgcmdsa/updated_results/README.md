# Updated Local Results

This directory is a compact archive for the local CMDSA reproduction results.
It stores only plotting inputs and the generated PDF report; dense BWR operator
matrices and raw log files are intentionally left as ignored generated data.

Files:

- `c5g7_local_runs.json`: local C5G7 timing and sweep rows.
- `bwr_spectral_diagnostics_with_fission.json`: BWR current-vs-paper metrics.
- `bwr_spectral_eigenvalues_with_fission.npz`: BWR eigenvalue arrays used for
  spectral plots.
- `updated_results.pdf`: standalone report generated from the archived data.

Refresh from local generated artifacts:

```bash
python3 ../run_all.py archive-results
python3 ../run_all.py updated-results
```
