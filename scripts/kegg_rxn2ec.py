#!/usr/bin/env python3
"""
Build rxn2ec.tsv by querying KEGG REST for KEGG Reaction IDs found in an SBML.

Requirements: network access. Use sparingly and cache results if you run often.

Usage:
  python scripts/kegg_rxn2ec.py --model models/drafts/creole_kbase.xml --out data/mappings/rxn2ec.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import time
import urllib.request
from pathlib import Path
from typing import Dict, List, Set, Tuple

from cobra.io import read_sbml_model


KEGG_RXN_RE = re.compile(r"R\d{5}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="KEGG REST -> rxn2ec.tsv")
    p.add_argument("--model", required=True, help="Input SBML path")
    p.add_argument("--out", required=True, help="Output TSV path")
    p.add_argument("--sleep", type=float, default=0.2, help="Sleep between requests (s)")
    return p.parse_args()


def extract_kegg_rids(model) -> Dict[str, Set[str]]:
    rxn_to_rids: Dict[str, Set[str]] = {}
    for rxn in model.reactions:
        ids: Set[str] = set()
        # Look in name and annotations
        if rxn.name:
            ids.update(KEGG_RXN_RE.findall(rxn.name))
        try:
            ann = rxn.annotation
        except Exception:
            ann = None
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
        blob = "\n".join(texts)
        ids.update(KEGG_RXN_RE.findall(blob))
        if ids:
            rxn_to_rids[rxn.id] = ids
    return rxn_to_rids


def kegg_reaction_to_ecs(rid: str) -> List[str]:
    url = f"https://rest.kegg.jp/get/rn:{rid}"
    with urllib.request.urlopen(url, timeout=8) as resp:
        text = resp.read().decode("utf-8", errors="ignore")
    ecs: Set[str] = set()
    for line in text.splitlines():
        if line.startswith("ENZYME") or line.startswith("            "):
            # lines with EC entries appear under ENZYME field; continuation lines start with spaces
            parts = line.split()
            for tok in parts[1:] if parts and parts[0] == "ENZYME" else parts:
                # ECs look like 1.1.1.1 or may be separated by spaces
                if tok.count(".") == 3:
                    ecs.add(tok.strip())
    return sorted(ecs)


def main() -> None:
    args = parse_args()
    model = read_sbml_model(args.model)
    rxn_to_rids = extract_kegg_rids(model)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    rows: List[Tuple[str, str]] = []
    for rxn_id, rids in sorted(rxn_to_rids.items()):
        ecs: Set[str] = set()
        for rid in sorted(rids):
            try:
                ecs.update(kegg_reaction_to_ecs(rid))
                time.sleep(args.sleep)
            except Exception as e:
                print(f"WARN: {rid}: {e}", file=sys.stderr)
        rows.append((rxn_id, ";".join(sorted(ecs))))

    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reaction_id", "ec"])
        w.writerows(rows)
    print(f"Wrote {len(rows)} reactions -> {out}")


if __name__ == "__main__":
    main()

