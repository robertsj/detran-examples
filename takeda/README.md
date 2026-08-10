# Takeda Examples

This folder contains Takeda Model 1 three-dimensional two-group benchmark
drivers from NEACRP-L-330, "3-D Neutron Transport Benchmarks."  The benchmark
has inserted- and withdrawn-control eigenvalue cases with reflective west,
south, and bottom boundaries and vacuum east, north, and top boundaries.

The historical scripts run the original SN/DD setup:

```sh
python3 takeda_inserted.py
python3 takeda_withdrawn.py
```

Native MOC runs are available through the shared Model 1 setup:

```sh
python3 run_moc.py --case inserted --equation scmoc
python3 run_moc.py --case inserted --equation ddmoc
python3 run_moc.py --case withdrawn --equation scmoc
```

`run_moc.py` defaults to a coarsened mesh with `--mesh-divisor 10`, so the
commands above are local smoke runs.  Use `--reference-mesh` or
`--mesh-divisor 1` for the original fine mesh from the legacy SN scripts.
`test_takeda.py` still includes the full legacy SN reference cases and can take
several minutes.

Useful options:

- `--equation scmoc|ddmoc` selects step-characteristics or diamond-difference
  MOC.
- `--polar`, `--azimuth`, and `--quadrature` select the angular product rule.
- `--spacing` controls the maximum average native track spacing.
- `--json` emits a machine-readable summary with keff, reference delta, mesh
  counts, sweep count, and scalar-flux statistics.

Source map:

| File | Purpose |
| --- | --- |
| `model.py` | Shared Takeda Model 1 mesh, material, reference, and solver builders. |
| `takeda_inserted.py` | Legacy SN/DD inserted-control wrapper. |
| `takeda_withdrawn.py` | Legacy SN/DD withdrawn-control wrapper. |
| `run_moc.py` | Native 3-D `scmoc`/`ddmoc` eigenvalue driver. |
| `test_takeda.py` | Regression checks for the legacy SN references and MOC smoke path. |
