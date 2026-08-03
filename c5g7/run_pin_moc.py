# pyexamples/c5g7/run_pin_moc.py
#
# Solves a reflected C5G7 pin eigenproblem with the native 2-D
# step-characteristics MOC sweeper.

try:
  from .moc_driver import parse_moc_args, print_summary, run_from_args, run_moc_case
except ImportError:
  from moc_driver import parse_moc_args, print_summary, run_from_args, run_moc_case


def run(**options):
  return run_moc_case("pin", **options)


def parse_args():
  return parse_moc_args(
    "Run native 2-D MOC on a C5G7 pin smoke case.",
    default_index=0,
  )


if __name__ == "__main__":
  print_summary(run_from_args(parse_args(), case="pin"))
