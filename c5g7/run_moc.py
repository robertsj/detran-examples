try:
  from .moc_driver import parse_moc_args, print_summary, run_from_args
except ImportError:
  from moc_driver import parse_moc_args, print_summary, run_from_args


def parse_args():
  return parse_moc_args(
    "Run native 2-D MOC on C5G7 pin, assembly, or core cases.",
    include_case=True,
    default_case="assembly",
  )


def run(args=None):
  if args is None:
    args = parse_args()
  return run_from_args(args)


if __name__ == "__main__":
  print_summary(run(parse_args()))
