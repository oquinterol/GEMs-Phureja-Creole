#!/usr/bin/env python3
"""
Extract rxn2ec.tsv (reaction_id\tec) from SBML annotations when ECs are present.

Looks for EC numbers in reaction annotations/notes (CVTerms or free text),
otherwise writes the reaction IDs with empty EC for manual completion.

Usage:
  python scripts/sbml_to_rxn2ec.py --model models/drafts/creole_kbase.xml --out data/mappings/rxn2ec.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Set

from cobra.io import read_sbml_model


EC_RE = re.compile(r"\b(\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-)\b")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="SBML -> rxn2ec.tsv extractor")
    p.add_argument("--model", required=True, help="Input SBML path")
    p.add_argument("--out", required=True, help="Output TSV path")
    return p.parse_args()


def extract_from_reaction(rxn) -> List[str]:
    ecs: Set[str] = set()
    # Try reaction annotation string if available
    ann = None
    try:
        ann = rxn.annotation  # dict-like in COBRApy
    except Exception:
        ann = None
    # Search in flattened string of annotations
    texts: List[str] = []
    if ann:
        def flatten(x):
            if isinstance(x, dict):
                for v in x.values():
                    yield from flatten(v)
            elif isinstance(x, (list, tuple, set)):
                for v in x:
                    yield from flatten(v)
            else:
                yield str(x)
        texts.extend(list(flatten(ann)))
    # Also check notes if present
    try:
        if rxn.notes:
            texts.append(str(rxn.notes))
    except Exception:
        pass
    blob = "\n".join(texts)
    for m in EC_RE.finditer(blob):
        ecs.add(m.group(0))
    return sorted(ecs)


def main() -> None:
    args = parse_args()
    model = read_sbml_model(args.model)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for rxn in model.reactions:
        ecs = extract_from_reaction(rxn)
        rows.append((rxn.id, ";".join(ecs)))
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reaction_id", "ec"])
        w.writerows(rows)
    print(f"Wrote {len(rows)} reactions -> {out}")


if __name__ == "__main__":
    main()

