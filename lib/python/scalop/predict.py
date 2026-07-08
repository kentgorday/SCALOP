# Assign canonical forms based on sequence-based length-independent clustering results
from __future__ import annotations

import argparse
import datetime
import json
import multiprocessing
import os
import pickle
import re
import sys
import time
from types import SimpleNamespace
from typing import Any

import numpy as np
from anarci import run_anarci

from ._types import Assignment
from .LoopGraft import export_structure
from .utils import SimpleFastaParser, getnumberedCDRloop, print_result, resns, write_txt

# Populated by ip() from the pickled cluster/PSSM database.
outm: dict[str, Any] = {}
clustercenters: dict[str, Any] = {}


def ip(scheme: str, definition: str, dbv: str) -> None:
    global outm
    global clustercenters
    scalop_path = os.path.split(__file__)[0]  # from scalop database

    tlist = []
    if dbv == "latest":
        for n in os.listdir(os.path.join(scalop_path, "database")):
            m = re.match(rf"{scheme}_{definition}_v(\d+)-(\d+).pickle", n)
            if m is None:
                continue
            nt = datetime.datetime.strptime("".join(m.groups()), "%Y%m")
            tlist.append((nt, n))
        if tlist == []:
            sys.stderr.write("Database is missing. Aborting!\n")
            sys.exit(1)

        fname = tlist[tlist.index(max(tlist))][1]
    else:
        fname = f"{scheme}_{definition}_v{dbv}.pickle"
    if not os.path.exists(os.path.join(scalop_path, "database", fname)):
        sys.stderr.write(f"Database {dbv} not found. Please review the database directory.")
        sys.exit(1)
    with open(os.path.join(scalop_path, "database", fname), "rb") as f:
        outm, clustercenters, _ = pickle.load(f, encoding="latin1")


def _score_licanonicalliPSSM(sequence: list[Any], pssm: dict[Any, Any]) -> float:
    score = float(0)
    for pos, res in sequence:
        if pos not in pssm or pssm[pos][resns.index(res.upper())] == np.nan:
            continue

        score += pssm[pos][resns.index(res)]
    return float(score) / float(len(sequence))


assignmentthresholds: dict[str, float] = {
    "H1": 0.5,
    "H2": -0.5,
    "L1": -0.5,
    "L2": -1,
    "L3": -1,
}  # legacy from 2017 PSSM
cdrs: dict[str, list[str]] = {"H": ["H1", "H2"], "L": ["L1", "L2", "L3"]}


def _assign(sequence: list[Any], cdr: str) -> list[str]:
    assignmentthreshold = assignmentthresholds[cdr]
    loopseq = "".join([x[1] for x in sequence]).upper()
    cdr = cdr.upper()
    if any([n not in resns for n in loopseq]):
        return [cdr, "None", "None", "None"]
    seqlen = len(sequence)
    _clusters = [
        cluster for cluster in sorted(outm[cdr]) if seqlen in outm[cdr][cluster]["Lengths"]
    ]  # clusters that have the sequence length
    # assign all length-8 L2 sequences (North et al. definition; or loops of the
    # majority length(s)) to 1 cluster
    if cdr == "L2":
        charac_seqlen = next(iter(outm["L2"].values()))["Lengths"]
        if seqlen in charac_seqlen:
            ID = next(iter(outm["L2"].keys()))
            center = clustercenters[cdr][ID]
            return [cdr, loopseq, ID, center]
        else:
            return [cdr, loopseq, "None", "None"]
    scores = [
        _score_licanonicalliPSSM(sequence, outm[cdr][cluster]["PSSM"])
        for cluster in sorted(outm[cdr])
        if seqlen in outm[cdr][cluster]["Lengths"]
    ]

    if any(np.greater(np.array(scores), assignmentthreshold)):
        pos = scores.index(max(scores))
        ID = _clusters[pos]
        center = clustercenters[cdr][ID]
        return [cdr, loopseq, ID, center]
    else:  # all were 0
        return [cdr, loopseq, "None", "None"]


def assignweb(args: argparse.Namespace) -> tuple[list[Assignment], str, str | bool]:
    seqs: list[tuple[str, str]] = []
    permitbuildstruct = False
    if os.path.isfile(args.seq):
        with open(args.seq) as f:
            _seqs = list(SimpleFastaParser(f))
        for seqid, seq in _seqs:
            if "/" in seq:
                seqh, seql = seq.split("/")
                seqh = re.sub(r"[\W]+", "", seqh)
                seql = re.sub(r"[\W]+", "", seql)
                seqs.append((seqid + "_1", seqh))
                seqs.append((seqid + "_2", seql))
            else:
                seqs.append((seqid, re.sub(r"[\W\/]+", "", seq)))
    else:  # Single sequence input
        permitbuildstruct = True
        if "/" in args.seq:  # must be Heavy-chain/Light-chain
            seqh, seql = args.seq.split("/")
            seqs = [("input_1", seqh), ("input_2", seql)]
        else:
            seqs = [("input", args.seq)]

    assignresults: list[Assignment] = [{} for _ in range(len(seqs))]

    # Number the heavy / light chain using ANARCI
    ncpu = multiprocessing.cpu_count()
    numberedseqs = run_anarci(seqs, scheme=args.scheme, ncpu=ncpu, assign_germline=False)

    # Import the right version of pssm database
    ip(args.scheme, args.definition, args.dbv)

    for i in range(len(numberedseqs[0])):
        seqname = numberedseqs[0][i][0]
        assert seqname == seqs[i][0]
        assignresults[i] = {"seqname": seqname, "input": seqs[i], "outputs": {}}
        # Find if ANARCI can number this chain (whether or not it is an antibody chain)
        if numberedseqs[2][i] is None:
            continue
        # Find whether heavy or light chain
        chain = "H" if numberedseqs[2][i][0]["chain_type"] == "H" else "L"

        for cdr in cdrs[chain]:
            # Frame the CDR region
            loop, _ = getnumberedCDRloop(
                numberedseqs[1][i][0][0], cdr, args.scheme, args.definition
            )

            # Assign canonical form
            assignresults[i]["outputs"].update({cdr: _assign(loop, cdr)})
    opf = write_txt(assignresults, args.outputdir, args.scheme, args.definition)
    outpdbf: str | bool = False

    if permitbuildstruct and args.structuref != "":
        _outpdbf = opf.replace(".txt", ".pdb")
        graftedstructs, is_model_generated = export_structure(args, assignresults, _outpdbf)
        if is_model_generated == 1:
            opf = write_txt(
                assignresults, args.outputdir, args.scheme, args.definition, graftedstructs
            )
            outpdbf = _outpdbf
    return assignresults, opf, outpdbf


def assigncmd(args: argparse.Namespace) -> list[Assignment]:
    seqs: list[tuple[str, str]] = []
    permitbuildstruct = False
    if os.path.isfile(args.seq):
        with open(args.seq) as f:
            _seqs = list(SimpleFastaParser(f))
        for seqid, seq in _seqs:
            if "/" in seq:
                seqh, seql = seq.split("/")
                seqh = re.sub(r"[\W]+", "", seqh)
                seql = re.sub(r"[\W]+", "", seql)
                seqs.append((seqid + "_1", seqh))
                seqs.append((seqid + "_2", seql))
            else:
                seqs.append((seqid, re.sub(r"[\W\/]+", "", seq)))
    else:  # Single sequence input
        permitbuildstruct = True
        if "/" in args.seq:  # must be Heavy-chain/Light-chain
            seqh, seql = args.seq.split("/")
            seqs = [("input_1", seqh), ("input_2", seql)]
        else:
            seqs = [("input", args.seq)]

    assignresults: list[Assignment] = [{} for _ in range(len(seqs))]

    # Number the heavy / light chain using ANARCI
    ncpu = multiprocessing.cpu_count()
    numberedseqs = run_anarci(seqs, scheme=args.scheme, ncpu=ncpu, assign_germline=False)

    # Import the right version of pssm database
    ip(args.scheme, args.definition, args.dbv)

    for i in range(len(numberedseqs[0])):
        seqname = numberedseqs[0][i][0]
        assert seqname == seqs[i][0]
        assignresults[i] = {"seqname": seqname, "input": seqs[i], "outputs": {}}
        # Find if ANARCI can number this chain (whether or not it is an antibody chain)
        if numberedseqs[2][i] is None:
            continue
        # Find whether heavy or light chain
        chain = "H" if numberedseqs[2][i][0]["chain_type"] == "H" else "L"

        for cdr in cdrs[chain]:
            # Frame the CDR region
            loop, _ = getnumberedCDRloop(
                numberedseqs[1][i][0][0], cdr, args.scheme, args.definition
            )

            # Assign canonical form
            assignresults[i]["outputs"].update({cdr: _assign(loop, cdr)})

    # Output format
    if not args.outputdir:
        print_result(assignresults, args.scheme, args.definition)
    elif args.outputdir == "":
        pass  # No output
    elif args.outputformat in ["txt", "csv"]:
        opf = write_txt(assignresults, args.outputdir, args.scheme, args.definition)
        if permitbuildstruct and args.structuref != "":
            outpdbf = opf.replace(".txt", ".pdb")
            graftedstructs, is_model_generated = export_structure(args, assignresults, outpdbf)
            if is_model_generated:
                opf = write_txt(
                    assignresults, args.outputdir, args.scheme, args.definition, graftedstructs
                )
        print(f"Write results to {opf}")
    elif args.outputformat == "json":
        opf = os.path.join(args.outputdir, f"{time.time()}.{args.outputformat}")
        if not os.path.exists(args.outputdir):
            os.mkdir(args.outputdir)
        with open(opf, "w") as f:
            json.dump(assignresults, f)
        print(f"Write results to {opf}")
    return assignresults


def assign(
    seq: str | list[str] | dict[str, str],
    scheme: str = "imgt",
    definition: str = "north",
    dbv: str = "latest",
    structuref: str = "",
    loopdb: str = "",
    hc: str = "",
    lc: str = "",
    blacklist: list[str] | None = None,
) -> list[Assignment]:
    """
    args: <sequence(s)> <numbering scheme> <cdr definition> <db version>
          <structure file> <loop database directory> <heavy chain ID>
          <light chain ID> <blacklisted PDB_CHAIN IDs>
    """

    if blacklist is None:
        blacklist = []
    seqs: list[tuple[str, str]] = []
    permitbuildstruct = False
    if isinstance(seq, str) and os.path.isfile(seq):
        with open(seq) as f:
            _seqs = list(SimpleFastaParser(f))

        for seqid, seq in _seqs:
            if "/" in seq:
                seqh, seql = seq.split("/")
                seqh = re.sub(r"[\W]+", "", seqh)
                seql = re.sub(r"[\W]+", "", seql)
                seqs.append((f"{seqid}_1", seqh))
                seqs.append((f"{seqid}_2", seql))
            else:
                seqs.append((seqid, re.sub(r"[\W\/]+", "", seq)))
    elif isinstance(seq, dict):  # in dict
        for seqn in seq:
            if "/" in seq[seqn]:  # must be Heavy-chain/Light-chain
                seqh, seql = seq[seqn].split("/")
                seqs.append((f"{seqn}_1", seqh))
                seqs.append((f"{seqn}_2", seql))
            else:
                seqs.append((seqn, seq[seqn]))

    elif isinstance(seq, list):  # in unnamed list
        for seqi in range(len(seq)):
            if "/" in seq[seqi]:  # must be Heavy-chain/Light-chain
                seqh, seql = seq[seqi].split("/")
                seqs.append((f"{seqi}_1", seqh))
                seqs.append((f"{seqi}_2", seql))
            else:
                seqs.append((str(seqi), seq[seqi]))

    else:  # Single sequence input
        permitbuildstruct = True
        if "/" in seq:  # must be Heavy-chain/Light-chain
            seqh, seql = seq.split("/")
            seqs = [("input_1", seqh), ("input_2", seql)]
        else:
            seqs = [("input", seq)]
    assignresults: list[Assignment] = [{} for _ in range(len(seqs))]
    # Number the heavy / light chain using ANARCI

    ncpu = multiprocessing.cpu_count()

    numberedseqs = run_anarci(seqs, scheme=scheme, ncpu=ncpu, assign_germline=False)

    # Import the right version of pssm database
    ip(scheme, definition, dbv)

    for i in range(len(numberedseqs[0])):
        # Find if ANARCI can number this chain (whether or not it is an antibody chain)
        if numberedseqs[2][i] is None:
            continue
        # Find whether heavy or light chain
        chain = "H" if numberedseqs[2][i][0]["chain_type"] == "H" else "L"
        seqname = numberedseqs[0][i][0]
        assert seqname == seqs[i][0]
        assignresults[i] = {"seqname": seqname, "input": seqs[i], "outputs": {}}
        for cdr in cdrs[chain]:
            # Frame the CDR region
            loop, _ = getnumberedCDRloop(numberedseqs[1][i][0][0], cdr, scheme, definition)

            # Assign canonical form
            assignresults[i]["outputs"].update({cdr: _assign(loop, cdr)})
    if permitbuildstruct and structuref != "":
        args = SimpleNamespace(
            scheme=scheme,
            definition=definition,
            dbv=dbv,
            structuref=os.path.abspath(structuref),
            loopdb=loopdb,
            hc=hc,
            lc=lc,
            blacklist=blacklist,
        )
        graftedstructs, _is_model_generated = export_structure(args, assignresults, "")

        for i in range(len(assignresults)):
            if "outputs" not in assignresults[i]:
                continue
            for cdr in assignresults[i]["outputs"]:
                if cdr in graftedstructs:
                    assignresults[i]["outputs"][cdr] += [graftedstructs[cdr]]
    return assignresults
