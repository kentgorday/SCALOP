"""Golden regression tests against SCALOP webserver output.

The reference files in tests/golden/webserver/ were downloaded from the OPIG
SCALOP webserver for the two example FASTAs across all six numbering/definition
combos. A one-combo cross-check confirmed the webserver's database == our bundled
`dbv=2019-10`, so these files double as a *verified-correct* oracle, not just a
snapshot of current behavior.

Two tiers:
  * Tier A  -- loop Sequence equality (DB-independent) plus oracle-free
    invariants (loop is a substring of the chain; the length encoded in the
    canonical id matches the loop length; a chain's CDRs are all-heavy or
    all-light). Breaks if numbering / region logic regresses.
  * Tier B  -- Canonical and Median equality. Breaks if database loading /
    scoring regresses.

Run with:

    pixi run test
"""

from pathlib import Path
from typing import Any

import pytest
from scalop._types import Assignment
from scalop.predict import assign

# A parsed result cell (loop_sequence, canonical_id, median_id) and the
# {seqname: {cdr: cell}} table shape used on both the expected and actual sides.
Cell = tuple[str, str, str]
Table = dict[str, dict[str, Cell]]

REPO = Path(__file__).resolve().parent.parent
GOLDEN_DIR = REPO / "tests" / "golden" / "webserver"
DBV = "2019-10"

# Golden-file name prefix -> the example FASTA the webserver was run on.
FASTA_FOR = {
    "sample": REPO / "example" / "sample_sequences.fasta",
    "kabat": REPO / "example" / "kabattrimmed.fasta",
}

HEAVY_CDRS = {"H1", "H2"}
LIGHT_CDRS = {"L1", "L2", "L3"}

GOLDEN_FILES = sorted(GOLDEN_DIR.glob("*.txt"))


def parse_webserver(path: Path) -> tuple[str, str, Table]:
    """Return (scheme, definition, {seqname: {cdr: (seq, canonical, median)}})."""
    lines = path.read_text().replace("\r", "").split("\n")
    scheme = lines[0].split(":", 1)[1].strip()
    definition = lines[1].split(":", 1)[1].strip()
    # lines[2] is the column header (Input/CDR/Sequence/Canonical/Median).
    table: Table = {}
    for line in lines[3:]:
        if not line.strip():
            continue
        seqname, cdr, seq, canonical, median = line.split("\t")
        table.setdefault(seqname, {})[cdr] = (seq, canonical, median)
    return scheme, definition, table


def actual_table(results: list[Assignment]) -> Table:
    """assign() results -> {seqname: {cdr: (seq, canonical, median)}} (numbered chains)."""
    table: Table = {}
    for r in results:
        for cdr, val in (r.get("outputs") or {}).items():
            _, seq, canonical, median = val[:4]
            table.setdefault(r["seqname"], {})[cdr] = (seq, canonical, median)
    return table


def _cells(table: Table) -> dict[tuple[str, str], Cell]:
    """Flatten to {(seqname, cdr): (seq, canonical, median)} for easy diffing."""
    return {(s, c): v for s, cdrs in table.items() for c, v in cdrs.items()}


@pytest.fixture(scope="module", params=GOLDEN_FILES, ids=lambda p: p.stem)
def combo(request: pytest.FixtureRequest) -> dict[str, Any]:
    path = request.param
    prefix = path.stem.split("_")[0]  # "sample" | "kabat"
    scheme, definition, expected = parse_webserver(path)
    results = assign(str(FASTA_FOR[prefix]), scheme=scheme, definition=definition, dbv=DBV)
    return {
        "name": path.stem,
        "expected": expected,
        "actual": actual_table(results),
        "results": results,
    }


def test_tier_a_sequences_and_invariants(combo: dict[str, Any]) -> None:
    expected, actual, results = combo["expected"], combo["actual"], combo["results"]

    # Same set of chains numbered as the webserver.
    only_ours = sorted(set(actual) - set(expected))
    only_web = sorted(set(expected) - set(actual))
    assert set(actual) == set(expected), (
        f"numbered-chain set differs -- only ours: {only_ours} ; only webserver: {only_web}"
    )

    exp, act = _cells(expected), _cells(actual)
    assert set(act) == set(exp), f"CDR key set differs: {set(act) ^ set(exp)}"

    # Tier A oracle: loop sequences match (DB-independent).
    seq_mismatch = {k: (act[k][0], exp[k][0]) for k in exp if act[k][0] != exp[k][0]}
    assert not seq_mismatch, f"loop sequence mismatches {{key: (ours, web)}}: {seq_mismatch}"

    # Oracle-free invariants.
    for r in results:
        outs = r.get("outputs") or {}
        if not outs:
            continue
        cdrs = set(outs)
        assert cdrs <= HEAVY_CDRS or cdrs <= LIGHT_CDRS, (
            f"{r['seqname']} has mixed heavy/light CDRs: {cdrs}"
        )
        chain_seq = r["input"][1]
        for cdr, val in outs.items():
            seq, canonical = val[1], val[2]
            if seq in ("", "None"):
                continue
            assert seq in chain_seq, (
                f"{r['seqname']} {cdr} loop {seq!r} is not a substring of the input chain"
            )
            if canonical != "None":
                lengths = {int(x) for x in canonical.split("-")[1].split(",")}
                assert len(seq) in lengths, (
                    f"{r['seqname']} {cdr} loop length {len(seq)} "
                    f"not among canonical {canonical} lengths {sorted(lengths)}"
                )


def test_tier_b_canonical_and_median(combo: dict[str, Any]) -> None:
    exp, act = _cells(combo["expected"]), _cells(combo["actual"])
    missing = ("<missing>", "<missing>", "<missing>")
    mism = {
        k: {"ours": act.get(k, missing)[1:], "web": exp[k][1:]}
        for k in exp
        if act.get(k, missing)[1:] != exp[k][1:]
    }
    assert not mism, f"canonical/median mismatches: {mism}"
