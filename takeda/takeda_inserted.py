# Takeda benchmark, Model 1, inserted control (case 1).
# Reference keff ~ 0.9624.

try:
  from .model import run_legacy
except ImportError:
  from model import run_legacy


def run(write_outputs=True):
  return run_legacy(
    "inserted",
    output_prefix="takeda_inserted",
    write_outputs=write_outputs,
  )["keff"]


if __name__ == "__main__":
  run()
