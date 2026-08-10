# Kobayashi 3-D MOC Examples

This folder contains native 3-D MOC fixed-source smoke drivers for the
Kobayashi void benchmark family.  The cases are one-group shielding problems
with a source block, near-void streaming regions, shield material, reflective
planes at low `x/y/z`, and vacuum planes at high `x/y/z`.

Run from this directory with the current pydetran install on `PYTHONPATH`:

```bash
export PYTHONPATH=/home/robertsj/Research/detran/install/python:$PYTHONPATH
export LD_LIBRARY_PATH=/home/robertsj/Research/detran/install/lib:$LD_LIBRARY_PATH
python3 run_moc.py --case 2 --equation scmoc
python3 run_moc.py --case 2 --equation ddmoc
```

The default mesh keeps the run short.  The Gnat reference inputs use one uniform
refinement of the listed coarse divisions; use `--refine 1` for that mesh.  The
driver currently reports solver size, material/source counts, sweep count, and
basic scalar-flux statistics.  The benchmark line samples are encoded in the
case metadata but are not yet interpolated into published-style line plots.

## Cases

| Case | Description |
| --- | --- |
| `1` | Volume source in a cubic shield. |
| `2` | Straight duct streaming problem. |
| `3` | Dog-legged duct streaming problem. |

## Provenance

The benchmark family is from Kobayashi et al., “3-D Radiation Transport
Benchmark Problems and Results for Simple Geometries with Void Region,”
OECD/NEA NSC/DOC(2000)4.  The material constants and coarse material maps here
mirror the public Gnat tutorial inputs for `kobayashi_1.i`, `kobayashi_2.i`, and
`kobayashi_3.i`.
