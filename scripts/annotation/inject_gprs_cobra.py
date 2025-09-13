#!/usr/bin/env python3
"""
Inject GPRs into an SBML model using simple EC-based mappings via COBRApy.

Inputs
- SBML draft (Path): e.g., models/drafts/creole_kbase.xml
- gene2ec.tsv: tab-separated with header: gene_id	ec
  - Multiple ECs per gene allowed, separated by ';' or ','
- rxn2ec.tsv: tab-separated with header: reaction_id	ec
  - Multiple ECs per reaction allowed, separated by ';' or ','

Logic (approximation)
- For each reaction, find all genes whose EC set intersects reaction EC set.
- Build GPR as OR over matching genes (isoenzymes). No complex subunit AND is inferred.
- Existing GPRs in the model are preserved if non-empty unless --overwrite is given.

Output
- Writes a new SBML with FBC gene associations at the given output path.

Usage
  python scripts/inject_gprs_cobra.py \
    --model models/drafts/creole_kbase.xml \
    --gene2ec data/mappings/gene2ec.tsv \
    --rxn2ec data/mappings/rxn2ec.tsv \
    --out models/curated/creole_with_gprs.xml

Requirements
- cobra (COBRApy)
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Set, List

import cobra
from cobra.io import read_sbml_model, write_sbml_model
from cobra.core.gene import Gene


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Inject EC-based GPRs into SBML via COBRApy")
    p.add_argument("--model", required=True, help="Input SBML path")
    p.add_argument("--gene2ec", required=True, help="TSV: gene_id\tec")
    p.add_argument("--rxn2ec", required=True, help="TSV: reaction_id\tec")
    p.add_argument("--out", required=True, help="Output SBML path")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing non-empty GPRs")
    p.add_argument("--gpr_table", default=None, help="Optional TSV (reaction_id\\tgpr) to override or supply curated rules")
    p.add_argument("--out_gpr_table", default=None, help="Optional path to write the applied GPR table (reaction_id\\tgpr)")
    return p.parse_args()


def _split_multi(value: str) -> List[str]:
    if value is None:
        return []
    # Accept ';' or ',' as separators, strip spaces
    parts = [v.strip() for v in value.replace(";", ",").split(",")]
    # Normalize ECs: strip optional 'EC:' prefix
    norm: List[str] = []
    for p in parts:
        if not p:
            continue
        if p.upper().startswith("EC:"):
            p = p.split(":", 1)[1]
        norm.append(p)
    return norm


def _rxn_aliases(rid: str) -> List[str]:
    """Return possible aliases for a reaction ID differing by optional 'R_' prefix.

    Many KBase exports include reaction IDs both as 'rxn00001_c0' and 'R_rxn00001_c0'.
    To be robust, we will consider both when looking up mappings.
    """
    rid = rid.strip()
    if rid.startswith("R_"):
        return [rid, rid[2:]]
    else:
        return [rid, f"R_{rid}"]


def load_gene2ec(path: Path) -> Dict[str, Set[str]]:
    mapping: Dict[str, Set[str]] = {}
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "gene_id" not in reader.fieldnames or "ec" not in reader.fieldnames:
            raise ValueError("gene2ec.tsv must have headers: gene_id\tec")
        for row in reader:
            gid = row["gene_id"].strip()
            ecs = set(_split_multi(row["ec"]))
            if not gid or not ecs:
                continue
            mapping.setdefault(gid, set()).update(ecs)
    return mapping


def load_rxn2ec(path: Path) -> Dict[str, Set[str]]:
    mapping: Dict[str, Set[str]] = {}
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "reaction_id" not in reader.fieldnames or "ec" not in reader.fieldnames:
            raise ValueError("rxn2ec.tsv must have headers: reaction_id\tec")
        for row in reader:
            rid = row["reaction_id"].strip()
            ecs = set(_split_multi(row["ec"]))
            if not rid or not ecs:
                continue
            for alias in _rxn_aliases(rid):
                mapping.setdefault(alias, set()).update(ecs)
    return mapping


def load_gpr_table(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "reaction_id" not in reader.fieldnames or "gpr" not in reader.fieldnames:
            raise ValueError("gpr_table must have headers: reaction_id\tgpr")
        for row in reader:
            rid = row["reaction_id"].strip()
            gpr = (row.get("gpr", "") or "").strip()
            if rid and gpr:
                # Store both alias forms
                for alias in _rxn_aliases(rid):
                    mapping[alias] = gpr
    return mapping


def main() -> None:
    args = parse_args()
    model_path = Path(args.model)
    out_path = Path(args.out)
    gene2ec = load_gene2ec(Path(args.gene2ec))
    rxn2ec = load_rxn2ec(Path(args.rxn2ec))
    overrides: Dict[str, str] = {}
    if args.gpr_table:
        overrides = load_gpr_table(Path(args.gpr_table))

    model: cobra.Model = read_sbml_model(str(model_path))

    # Build reverse index: EC -> genes
    ec2genes: Dict[str, Set[str]] = {}
    for g, ecs in gene2ec.items():
        for ec in ecs:
            ec2genes.setdefault(ec, set()).add(g)

    injected = 0
    preserved = 0
    skipped = 0
    missing = 0

    applied_rows: List[List[str]] = []

    for rxn in model.reactions:
        # Try direct and alias lookups for reaction ECs
        rxn_ecs: Set[str] | None = None
        for alias in _rxn_aliases(rxn.id):
            rxn_ecs = rxn2ec.get(alias)
            if rxn_ecs:
                break
        if not rxn_ecs and rxn.id not in overrides:
            missing += 1
            continue
        # Collect genes whose EC overlaps
        genes_for_rxn: Set[str] = set()
        if rxn_ecs:
            for ec in rxn_ecs:
                genes_for_rxn.update(ec2genes.get(ec, set()))
        if not genes_for_rxn:
            skipped += 1
            continue

        # Determine if we should overwrite existing rule
        if rxn.gene_reaction_rule and not args.overwrite and rxn.id not in overrides:
            preserved += 1
            continue

        # Determine rule: override takes precedence, else OR of genes
        if rxn.id in overrides:
            rule = overrides[rxn.id]
            # Ensure referenced genes exist; if missing, add to model.genes
            for tok in re.findall(r"[A-Za-z0-9_.:-]+", rule):
                if tok.lower() in {"and", "or"}:
                    continue
                if tok not in model.genes:
                    model.genes.add(Gene(tok))
        else:
            # Build OR rule
            for gid in genes_for_rxn:
                if gid not in model.genes:
                    model.genes.add(Gene(gid))
            rule = " or ".join(sorted(genes_for_rxn))
        rxn.gene_reaction_rule = rule
        applied_rows.append([rxn.id, rule])
        injected += 1

    write_sbml_model(model, str(out_path))

    print(
        f"Injected: {injected}, Preserved: {preserved}, Skipped(no genes): {skipped}, Missing(rxn2ec): {missing}",
        flush=True,
    )

    # Optionally write applied GPR table
    if args.out_gpr_table:
        out_gpr = Path(args.out_gpr_table)
        out_gpr.parent.mkdir(parents=True, exist_ok=True)
        with out_gpr.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["reaction_id", "gpr"])
            for rid, gpr in sorted(applied_rows):
                w.writerow([rid, gpr])


if __name__ == "__main__":
    main()
