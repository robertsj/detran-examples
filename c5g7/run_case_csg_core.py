"""Generate C5G7 CSG quarter-core benchmark artifacts."""

try:
  from .csg_c5g7 import main
except ImportError:
  from csg_c5g7 import main


if __name__ == "__main__":
  main(default_case="core", default_index=0, default_stem="csg_core_0")
