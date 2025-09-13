#!/usr/bin/env python3
"""
Build a preliminary rxn2gpr.tsv (reaction_id\tgpr) by intersecting rxn2ec.tsv and gene2ec.tsv.

This allows manual curation of multi-subunit (AND) rules before injecting into SBML.

Usage:
  python scripts/build_rxn2gpr_table.py \
    --gene2ec data/mappings/gene2ec.tsv \
    --rxn2ec data/mappings/rxn2ec.tsv \
    --out data/mappings/rxn2gpr.tsv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Set, List


def _split_multi(value: str) -> List[str]:
    if value is None:
        return []
    parts = [v.strip() for v in value.replace(";", ",").split(",")]
    return [p for p in parts if p]


def load_gene2ec(path: Path) -> Dict[str, Set[str]]:
    mapping: Dict[str, Set[str]] = {}
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = row.get("gene_id", "").strip()
            ecs = set(_split_multi(row.get("ec", "")))
            if gid and ecs:
                mapping.setdefault(gid, set()).update(ecs)
    return mapping


def load_rxn2ec(path: Path) -> Dict[str, Set[str]]:
    mapping: Dict[str, Set[str]] = {}
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rid = row.get("reaction_id", "").strip()
            ecs = set(_split_multi(row.get("ec", "")))
            if rid and ecs:
                mapping.setdefault(rid, set()).update(ecs)
    return mapping


def main() -> None:
    p = argparse.ArgumentParser(description="Build preliminary rxn2gpr table from EC mappings")
    p.add_argument("--gene2ec", required=True)
    p.add_argument("--rxn2ec", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    g2e = load_gene2ec(Path(args.gene2ec))
    r2e = load_rxn2ec(Path(args.rxn2ec))

    # Build EC -> genes map
    e2g: Dict[str, Set[str]] = {}
    for g, ecs in g2e.items():
        for ec in ecs:
            e2g.setdefault(ec, set()).add(g)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for rid, ecs in sorted(r2e.items()):
        genes: Set[str] = set()
        for ec in ecs:
            genes.update(e2g.get(ec, set()))
        gpr = " or ".join(sorted(genes)) if genes else ""
        rows.append((rid, gpr))

    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reaction_id", "gpr"])
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows -> {out_path}")


if __name__ == "__main__":
    main()

