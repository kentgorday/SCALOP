"""Shared type definitions for scalop."""

from __future__ import annotations

from typing import TypedDict


class Assignment(TypedDict, total=False):
    """One input chain's canonical-form assignment result.

    An entry is empty (``{}``) when ANARCI could not number the chain. Otherwise:

    * ``seqname`` -- the input identifier
    * ``input``   -- ``(name, sequence)``
    * ``outputs`` -- ``{CDR: [cdr, loop_seq, canonical_id, median_id]}`` (a fifth
      grafted-loop element is appended when structure grafting runs)
    """

    seqname: str
    input: tuple[str, str]
    outputs: dict[str, list[str]]
