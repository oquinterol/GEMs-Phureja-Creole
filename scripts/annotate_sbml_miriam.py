#!/usr/bin/env python3
"""
Annotate an SBML model with MIRIAM/Identifiers.org URIs (reactions and genes)
using local mapping files (eggNOG annotations and rxn2ec.tsv). Writes SBML L3+FBC.

Adds:
- Reactions: EC URIs from rxn2ec.tsv; KEGG reaction URIs found in reaction name/annotation.
- Genes (geneProducts): KO, GO, EC (from eggNOG .annotations), and KEGG reactions linked in eggNOG.

Notes
- This uses COBRApy's annotation dict (BQB_IS) which the SBML writer exports as bqbiol:is.
- No network calls are made; if KEGG lookups are required, use scripts/kegg_rxn2ec.py first.

Usage
  python scripts/annotate_sbml_miriam.py \
    --model models/curated/creole_with_gprs.xml \
    --emapper reports/emapper_creole.emapper.annotations \
    --rxn2ec data/mappings/rxn2ec.tsv \
    --out models/curated/creole_with_gprs_annot.xml
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Set

import cobra
from cobra.io import read_sbml_model, write_sbml_model
from xml.etree import ElementTree as ET


KEGG_RXN_RE = re.compile(r"R\d{5}")
EC_RE = re.compile(r"\b(\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-)\b")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Annotate SBML with MIRIAM URIs using eggNOG + rxn2ec")
    p.add_argument("--model", required=True, help="Input SBML path (with GPRs)")
    p.add_argument("--emapper", required=True, help="eggNOG .annotations file (TSV)")
    p.add_argument("--rxn2ec", required=True, help="TSV reaction_id\tec")
    p.add_argument("--out", required=True, help="Output SBML path")
    return p.parse_args()


def _split_multi(value: str) -> List[str]:
    if not value:
        return []
    return [v.strip() for v in re.split(r"[;,]", value) if v.strip()]


def load_rxn2ec(path: Path) -> Dict[str, Set[str]]:
    mapping: Dict[str, Set[str]] = {}
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        assert reader.fieldnames and "reaction_id" in reader.fieldnames and "ec" in reader.fieldnames, "rxn2ec header required"
        for row in reader:
            rid = row["reaction_id"].strip()
            ecs = {ec if not ec.upper().startswith("EC:") else ec.split(":", 1)[1] for ec in _split_multi(row["ec"]) if ec}
            if not rid or not ecs:
                continue
            # accept both aliases (with and without R_)
            keys = [rid, rid[2:]] if rid.startswith("R_") else [rid, f"R_{rid}"]
            for k in keys:
                mapping.setdefault(k, set()).update(ecs)
    return mapping


def parse_emapper_annotations(path: Path) -> Dict[str, Dict[str, Set[str]]]:
    """Return gene -> {'ko': set, 'go': set, 'ec': set, 'kegg_rxn': set} from emapper .annotations."""
    gmap: Dict[str, Dict[str, Set[str]]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        header: List[str] | None = None
        for line in fh:
            if line.startswith("#"):
                # header line typically starts with '#', next line is header or columns embedded
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                # Normalize column indices
                col_idx = {name: i for i, name in enumerate(header)}
                # Minimal columns expected
                # Handle different emapper schemas by probing known names
                cols = {
                    "query": col_idx.get("query_name", col_idx.get("#query_name", 0)),
                    "GOs": col_idx.get("GOs"),
                    "EC": col_idx.get("EC"),
                    "KO": col_idx.get("KEGG_ko"),
                    "KEGG_RXN": col_idx.get("KEGG_Reaction"),
                }
                # store for closure
                _cols = cols
                continue
            cols = _cols  # type: ignore[name-defined]
            if not parts or len(parts) <= cols["query"]:
                continue
            gid = parts[cols["query"]].strip()
            if not gid:
                continue
            entry = gmap.setdefault(gid, {"ko": set(), "go": set(), "ec": set(), "kegg_rxn": set()})
            # GO terms
            if cols["GOs"] is not None and cols["GOs"] < len(parts):
                for go in _split_multi(parts[cols["GOs"]]):
                    if go.startswith("GO:"):
                        entry["go"].add(go)
            # EC numbers
            if cols["EC"] is not None and cols["EC"] < len(parts):
                for ec in _split_multi(parts[cols["EC"]]):
                    if EC_RE.fullmatch(ec):
                        entry["ec"].add(ec)
            # KO terms
            if cols["KO"] is not None and cols["KO"] < len(parts):
                for ko in _split_multi(parts[cols["KO"]]):
                    if ko.startswith("ko:"):
                        ko = ko.split(":", 1)[1]
                    if ko and re.fullmatch(r"K\d{5}", ko):
                        entry["ko"].add(ko)
            # KEGG reactions linked to the gene (optional)
            if cols["KEGG_RXN"] is not None and cols["KEGG_RXN"] < len(parts):
                for rx in _split_multi(parts[cols["KEGG_RXN"]]):
                    for m in KEGG_RXN_RE.finditer(rx):
                        entry["kegg_rxn"].add(m.group(0))
    return gmap


def add_uri(ann: dict, uri: str) -> None:
    vals = ann.setdefault("BQB_IS", [])
    if uri not in vals:
        vals.append(uri)


def main() -> None:
    args = parse_args()
    model_path = Path(args.model)
    out_path = Path(args.out)

    model: cobra.Model = read_sbml_model(str(model_path))
    rxn2ec = load_rxn2ec(Path(args.rxn2ec))
    gene_ann = parse_emapper_annotations(Path(args.emapper))

    # Annotate reactions with EC and KEGG reaction IDs
    for rxn in model.reactions:
        ann = dict(rxn.annotation or {})
        ecs = set()
        # from rxn2ec table
        ecs |= rxn2ec.get(rxn.id, set())
        # parse existing annotations and name for KEGG RIDs
        blob = [rxn.name or ""]
        try:
            def flatten(x):
                if isinstance(x, dict):
                    for v in x.values():
                        yield from flatten(v)
                elif isinstance(x, (list, tuple, set)):
                    for v in x:
                        yield from flatten(v)
                else:
                    yield str(x)
            blob.extend(list(flatten(rxn.annotation)))
        except Exception:
            pass
        text = "\n".join(blob)
        for rid in KEGG_RXN_RE.findall(text):
            add_uri(ann, f"https://identifiers.org/kegg.reaction/{rid}")
        for ec in ecs:
            add_uri(ann, f"https://identifiers.org/ec-code/{ec}")
        # Assign SBO for typical pseudo-reactions if missing
        if not getattr(rxn, "sbo", None):
            if rxn.id.startswith("EX_"):
                rxn.sbo = "SBO:0000627"  # exchange
            elif rxn.id.startswith("DM_"):
                rxn.sbo = "SBO:0000628"  # demand
            elif "biomass" in rxn.id:
                rxn.sbo = "SBO:0000629"  # biomass production
        rxn.annotation = ann

    # Annotate genes with KO, GO, EC, and KEGG reactions
    for gene in model.genes:
        ann = dict(gene.annotation or {})
        ginfo = gene_ann.get(gene.id)
        if not ginfo:
            continue
        for ko in sorted(ginfo["ko"]):
            add_uri(ann, f"https://identifiers.org/kegg.orthology/{ko}")
        for go in sorted(ginfo["go"]):
            add_uri(ann, f"https://identifiers.org/go/{go.split(':',1)[1]}")
        for ec in sorted(ginfo["ec"]):
            add_uri(ann, f"https://identifiers.org/ec-code/{ec}")
        for rid in sorted(ginfo["kegg_rxn"]):
            add_uri(ann, f"https://identifiers.org/kegg.reaction/{rid}")
        gene.annotation = ann

    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_sbml_model(model, str(out_path))
    # Post-process: add GO to compartments via XML (COBRApy lacks compartment annotations API)
    add_go_to_compartments(out_path)
    print(f"Annotated model written to: {out_path}")


def add_go_to_compartments(sbml_path: Path) -> None:
    comp_go = {
        # Mapping derived from KBase documentation provided by user
        "c": "GO:0005737",  # Cytosol
        "d": "GO:0009570",  # Chloroplast stroma (Stroma)
        "e": "GO:0048046",  # Apoplast (Extracellular in plants)
        "g": "GO:0005794",  # Golgi apparatus
        "j": "GO:0005743",  # Mitochondrial inner membrane
        "m": "GO:0005739",  # Mitochondrion
        "n": "GO:0005634",  # Nucleus
        "r": "GO:0005783",  # Endoplasmic reticulum
        "v": "GO:0005773",  # Vacuole
        "w": "GO:0009505",  # Plant-type cell wall
        "x": "GO:0005777",  # Peroxisome
        # keep plastid mapping for potential 'p' compartment ids
        "p": "GO:0009536",  # Plastid
    }
    tree = ET.parse(str(sbml_path))
    root = tree.getroot()
    # namespaces (do not use register_namespace here to avoid modifying prefixes globally)
    RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    BQBIOL_NS = "http://biomodels.net/biology-qualifiers/"

    def q(ns, tag):
        return f"{{{ns}}}{tag}"

    # Find listOfCompartments
    for comp in root.findall('.//{http://www.sbml.org/sbml/level3/version1/core}compartment'):
        cid = comp.attrib.get("id")
        if not cid:
            continue
        prefix = cid[0]
        go_id = comp_go.get(prefix)
        if not go_id:
            continue
        uri = f"https://identifiers.org/{go_id}"

        # Check if already present
        already = False
        for li in comp.findall(f'.//{q(RDF_NS,"li")}'):
            if li.attrib.get(q(RDF_NS, "resource")) == uri or li.attrib.get("resource") == uri:
                already = True
                break
        if already:
            continue

        # Build annotation block (append or create)
        ann = comp.find('{http://www.sbml.org/sbml/level3/version1/core}annotation')
        if ann is None:
            ann = ET.SubElement(comp, '{http://www.sbml.org/sbml/level3/version1/core}annotation')
        rdf = ET.SubElement(ann, q(RDF_NS, 'RDF'))
        # set namespaces on this rdf element for readability (optional for many parsers)
        rdf.set("xmlns:rdf", RDF_NS)
        rdf.set("xmlns:bqbiol", BQBIOL_NS)
        desc = ET.SubElement(rdf, q(RDF_NS, 'Description'))
        desc.set(q(RDF_NS, 'about'), f"#{cid}")
        is_el = ET.SubElement(desc, q(BQBIOL_NS, 'is'))
        bag = ET.SubElement(is_el, q(RDF_NS, 'Bag'))
        li = ET.SubElement(bag, q(RDF_NS, 'li'))
        li.set(q(RDF_NS, 'resource'), uri)

    tree.write(str(sbml_path), encoding='utf-8', xml_declaration=True)


if __name__ == "__main__":
    main()
