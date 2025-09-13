#!/usr/bin/env python3
"""
Convert eggNOG-mapper output to gene2ec.tsv (gene_id\tec).

Accepts either .emapper.annotations (tab) or .emapper.csv.
Will look for ECs in known columns and/or parse free-text annotations.

Usage:
  python scripts/emapper_to_gene2ec.py \
    --input reports/emapper.emapper.annotations \
    --out data/mappings/gene2ec.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Dict, Set


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="eggNOG-mapper -> gene2ec.tsv")
    p.add_argument("--input", required=True, help="eggNOG-mapper .annotations or .csv")
    p.add_argument("--out", required=True, help="Output TSV path")
    return p.parse_args()


def extract_ec_tokens(text: str) -> List[str]:
    ecs: Set[str] = set()
    if not text:
        return []
    # Simple EC regex-less extraction: tokens like '1.1.1.1' (digits and dots, 4 fields)
    # Avoid false positives by requiring exactly 4 components and digits or '-' in each.
    for part in text.replace(";", ",").replace("|", ",").split(","):
        token = part.strip()
        if not token:
            continue
        comps = token.split(".")
        if len(comps) != 4:
            continue
        ok = True
        for c in comps:
            if c != '-' and not c.isdigit():
                ok = False
                break
        if ok:
            ecs.add(token)
    return sorted(ecs)


def from_annotations(path: Path) -> Dict[str, List[str]]:
    gene_ecs: Dict[str, Set[str]] = {}
    # eggNOG-mapper annotations are tab-delimited with hash-prefixed header lines
    # We read manually line by line
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if not cols:
                continue
            gene = cols[0]
            # Heuristic: ECs often appear in last columns or in the annotations col
            ecs: Set[str] = set()
            for c in cols[1:]:
                for ec in extract_ec_tokens(c):
                    ecs.add(ec)
            if ecs:
                gene_ecs.setdefault(gene, set()).update(ecs)
    return {g: sorted(list(v)) for g, v in gene_ecs.items()}


def from_csv(path: Path) -> Dict[str, List[str]]:
    gene_ecs: Dict[str, Set[str]] = {}
    with path.open("r", newline="", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        # Common column names (vary by version)
        possible_cols = [
            "EC", "EC_number", "ECs", "EC_Number(s)", "EC (eggnog)"
        ]
        ec_cols = [c for c in possible_cols if c in reader.fieldnames]
        for row in reader:
            gene = row.get("query", row.get("#query", row.get("gene", "")).strip())
            ecs: Set[str] = set()
            for cname in ec_cols:
                ecs.update(extract_ec_tokens(row.get(cname, "")))
            # Fallback: scan all fields
            if not ecs:
                for v in row.values():
                    ecs.update(extract_ec_tokens(v or ""))
            if gene and ecs:
                gene_ecs.setdefault(gene, set()).update(ecs)
    return {g: sorted(list(v)) for g, v in gene_ecs.items()}


def write_gene2ec(g2e: Dict[str, List[str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "ec"])
        for g, ecs in sorted(g2e.items()):
            w.writerow([g, ";".join(ecs)])


def main() -> None:
    args = parse_args()
    p = Path(args.input)
    if not p.exists():
        print(f"ERROR: input not found: {p}", file=sys.stderr)
        sys.exit(2)
    if p.suffix == ".csv":
        g2e = from_csv(p)
    else:
        g2e = from_annotations(p)
    write_gene2ec(g2e, Path(args.out))
    print(f"Wrote {len(g2e)} gene EC mappings -> {args.out}")


if __name__ == "__main__":
    main()

