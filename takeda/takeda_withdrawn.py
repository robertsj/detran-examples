# Takeda benchmark, Model 1, withdrawn control (case 2).
# Reference keff ~ 0.9780.

try:
  from .model import run_legacy
except ImportError:
  from model import run_legacy


def run(write_outputs=True):
  return run_legacy(
    "withdrawn",
    output_prefix="takeda_withdrawn",
    write_outputs=write_outputs,
  )["keff"]


if __name__ == "__main__":
  run()
