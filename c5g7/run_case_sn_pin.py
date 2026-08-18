"""Generate C5G7 SN/DD pin-cell benchmark artifacts."""

try:
  from .sn_driver import main
except ImportError:
  from sn_driver import main


if __name__ == "__main__":
  main(default_case="pin", default_index=0, default_stem="sn_pin_0_optimized")
