# pyexamples/c5g7/run_core_moc.py
#
# Solves the C5G7 core eigenproblem with the native 2-D
# step-characteristics MOC sweeper.

try:
  from .moc_driver import parse_moc_args, print_summary, run_from_args, run_moc_case
except ImportError:
  from moc_driver import parse_moc_args, print_summary, run_from_args, run_moc_case


def run(**options):
  return run_moc_case("core", **options)


def parse_args():
  return parse_moc_args(
    "Run native 2-D MOC on the C5G7 core case.",
    default_index=0,
  )


if __name__ == "__main__":
  print_summary(run_from_args(parse_args(), case="core"))
