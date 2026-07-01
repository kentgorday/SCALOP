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

The pristine source only runs under the period-appropriate env, so:

    pixi run -e legacy test
"""
from pathlib import Path

import pytest

from scalop.predict import assign

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


def parse_webserver(path):
    """Return (scheme, definition, {seqname: {cdr: (seq, canonical, median)}})."""
    lines = path.read_text().replace("\r", "").split("\n")
    scheme = lines[0].split(":", 1)[1].strip()
    definition = lines[1].split(":", 1)[1].strip()
    # lines[2] is the column header (Input/CDR/Sequence/Canonical/Median).
    table = {}
    for line in lines[3:]:
        if not line.strip():
            continue
        seqname, cdr, seq, canonical, median = line.split("\t")
        table.setdefault(seqname, {})[cdr] = (seq, canonical, median)
    return scheme, definition, table


def actual_table(results):
    """assign() results -> {seqname: {cdr: (seq, canonical, median)}} (numbered chains)."""
    table = {}
    for r in results:
        for cdr, val in (r.get("outputs") or {}).items():
            _, seq, canonical, median = val[:4]
            table.setdefault(r["seqname"], {})[cdr] = (seq, canonical, median)
    return table


def _cells(table):
    """Flatten to {(seqname, cdr): (seq, canonical, median)} for easy diffing."""
    return {(s, c): v for s, cdrs in table.items() for c, v in cdrs.items()}


@pytest.fixture(scope="module", params=GOLDEN_FILES, ids=lambda p: p.stem)
def combo(request):
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


def test_tier_a_sequences_and_invariants(combo):
    expected, actual, results = combo["expected"], combo["actual"], combo["results"]

    # Same set of chains numbered as the webserver.
    assert set(actual) == set(expected), (
        "numbered-chain set differs -- only ours: %s ; only webserver: %s"
        % (sorted(set(actual) - set(expected)), sorted(set(expected) - set(actual)))
    )

    exp, act = _cells(expected), _cells(actual)
    assert set(act) == set(exp), "CDR key set differs: %s" % (set(act) ^ set(exp))

    # Tier A oracle: loop sequences match (DB-independent).
    seq_mismatch = {k: (act[k][0], exp[k][0]) for k in exp if act[k][0] != exp[k][0]}
    assert not seq_mismatch, "loop sequence mismatches {key: (ours, web)}: %s" % seq_mismatch

    # Oracle-free invariants.
    for r in results:
        outs = r.get("outputs") or {}
        if not outs:
            continue
        cdrs = set(outs)
        assert cdrs <= HEAVY_CDRS or cdrs <= LIGHT_CDRS, (
            "%s has mixed heavy/light CDRs: %s" % (r["seqname"], cdrs)
        )
        chain_seq = r["input"][1]
        for cdr, val in outs.items():
            seq, canonical = val[1], val[2]
            if seq in ("", "None"):
                continue
            assert seq in chain_seq, (
                "%s %s loop %r is not a substring of the input chain" % (r["seqname"], cdr, seq)
            )
            if canonical != "None":
                lengths = {int(x) for x in canonical.split("-")[1].split(",")}
                assert len(seq) in lengths, (
                    "%s %s loop length %d not among canonical %s lengths %s"
                    % (r["seqname"], cdr, len(seq), canonical, sorted(lengths))
                )


def test_tier_b_canonical_and_median(combo):
    exp, act = _cells(combo["expected"]), _cells(combo["actual"])
    missing = ("<missing>", "<missing>", "<missing>")
    mism = {
        k: {"ours": act.get(k, missing)[1:], "web": exp[k][1:]}
        for k in exp
        if act.get(k, missing)[1:] != exp[k][1:]
    }
    assert not mism, "canonical/median mismatches: %s" % mism
