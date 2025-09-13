#!/usr/bin/env python3
"""
Extract rxn2ec.tsv (reaction_id\tec) directly from SBML XML without COBRA/libSBML.

This scans each <reaction> element and collects EC numbers found in its annotations
(e.g., rdf:li rdf:resource="https://identifiers.org/ec-code/1.2.3.4") or any text
within the reaction subtree that matches an EC pattern (d.d.d.d, allowing '-').

Usage:
  python scripts/sbml_to_rxn2ec_xml.py --model models/drafts/creole_kbase.xml --out data/mappings/rxn2ec.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from xml.etree import ElementTree as ET


EC_RE = re.compile(r"\b(\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-)\b")


def local(tag: str) -> str:
    return tag.split('}', 1)[-1] if '}' in tag else tag


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="SBML XML -> rxn2ec.tsv extractor (no dependencies)")
    p.add_argument("--model", required=True, help="Input SBML path")
    p.add_argument("--out", required=True, help="Output TSV path")
    return p.parse_args()


def extract_rxn2ec(xml_path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    # iterparse to keep memory low
    context = ET.iterparse(str(xml_path), events=("start", "end"))
    current_rxn_id: str | None = None
    collected_ecs: set[str] | None = None

    for event, elem in context:
        tag = local(elem.tag)
        if event == "start" and tag == "reaction":
            current_rxn_id = elem.attrib.get("id") or elem.attrib.get("metaid")
            collected_ecs = set()
        elif event == "end" and tag == "reaction":
            if current_rxn_id is not None and collected_ecs is not None:
                rows.append((current_rxn_id, ";".join(sorted(collected_ecs))))
            # clear the element to free memory
            elem.clear()
            current_rxn_id = None
            collected_ecs = None
        elif current_rxn_id is not None and collected_ecs is not None:
            # Within a reaction subtree: look for ECs in attributes and text
            # Check attributes (e.g., rdf:resource)
            for k, v in elem.attrib.items():
                if v:
                    # identifiers.org/ec-code/<EC>
                    if "identifiers.org/ec-code/" in v:
                        ec = v.rsplit("/", 1)[-1]
                        if EC_RE.fullmatch(ec):
                            collected_ecs.add(ec)
                    else:
                        for m in EC_RE.finditer(v):
                            collected_ecs.add(m.group(0))
            # Check text content
            if elem.text:
                for m in EC_RE.finditer(elem.text):
                    collected_ecs.add(m.group(0))
            if elem.tail:
                for m in EC_RE.finditer(elem.tail):
                    collected_ecs.add(m.group(0))
    return rows


def main() -> None:
    args = parse_args()
    in_path = Path(args.model)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows = extract_rxn2ec(in_path)
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reaction_id", "ec"])
        w.writerows(rows)
    print(f"Wrote {len(rows)} reactions -> {out_path}")


if __name__ == "__main__":
    main()

