"""Generate C5G7 SN/DD quarter-core benchmark artifacts."""

try:
  from .sn_driver import main
except ImportError:
  from sn_driver import main


if __name__ == "__main__":
  main(default_case="core", default_index=0, default_stem="sn_core_0_optimized")
