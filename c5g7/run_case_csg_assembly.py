"""Generate C5G7 CSG assembly benchmark artifacts."""

try:
  from .csg_c5g7 import main
except ImportError:
  from csg_c5g7 import main


if __name__ == "__main__":
  main(default_case="assembly", default_index=0, default_stem="csg_assembly_0")
