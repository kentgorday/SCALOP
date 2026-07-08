"""Smoke tests: run SCALOP on the bundled example inputs and assert only that it
produces output.

These intentionally check *that* outputs exist, not *what* they are. Their jobs:
  1. prove the environment can run the tool end-to-end (ANARCI -> hmmscan ->
     canonical-form assignment), and
  2. be the scaffold we later grow into value-checking golden/regression tests.

Run with:

    pixi run test
"""

from pathlib import Path

import pytest
from scalop.predict import assign

REPO = Path(__file__).resolve().parent.parent
EXAMPLES = REPO / "example"

# Pin the database snapshot so runs are reproducible (never 'latest').
DBV = "2019-10"

EXAMPLE_FASTAS = [
    EXAMPLES / "sample_sequences.fasta",  # clean mouse heavy/light chains
    EXAMPLES / "kabattrimmed.fasta",  # messy chicken chains (edge cases)
]


@pytest.mark.parametrize("fasta", EXAMPLE_FASTAS, ids=lambda p: p.name)
def test_assign_runs_on_example_fasta(fasta: Path) -> None:
    assert fasta.is_file(), f"missing example input: {fasta}"

    results = assign(str(fasta), dbv=DBV)

    # One result entry per chain, each a dict. Just assert we got something back.
    assert isinstance(results, list) and results, "assign() returned no results"
    # At least one chain should have produced CDR canonical-form outputs.
    assert any(r.get("outputs") for r in results), "no CDR outputs were produced"
