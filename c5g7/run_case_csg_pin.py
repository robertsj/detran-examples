"""Generate C5G7 CSG pin-cell benchmark artifacts."""

try:
  from .csg_c5g7 import main
except ImportError:
  from csg_c5g7 import main


if __name__ == "__main__":
  main(default_case="pin", default_index=0, default_stem="csg_pin_0")
