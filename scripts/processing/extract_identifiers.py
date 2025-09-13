#!/usr/bin/env python3
"""
Extract identifiers from SBML models for API mapping.

Creates lists of identifiers that can be used with api_database_mapper.py:
- KEGG compound IDs from metabolites
- KEGG reaction IDs from reactions  
- EC numbers from reactions
- Gene names for UniProt searches

Usage:
  python scripts/extract_identifiers.py --model models/current/creole_unified.xml --output-dir data/mappings/
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Set
import cobra
from cobra.io import read_sbml_model


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract identifiers from SBML model for API mapping")
    p.add_argument("--model", required=True, help="Input SBML model path")
    p.add_argument("--output-dir", required=True, help="Output directory for identifier lists")
    return p.parse_args()


def extract_kegg_compounds(model: cobra.Model) -> Set[str]:
    """Extract KEGG compound IDs from metabolite annotations."""
    kegg_compounds = set()
    
    for metabolite in model.metabolites:
        # Check annotations for KEGG compound IDs
        if hasattr(metabolite, 'annotation'):
            for key, value in metabolite.annotation.items():
                if 'kegg.compound' in key.lower():
                    if isinstance(value, list):
                        for v in value:
                            if isinstance(v, str) and v.startswith('C'):
                                kegg_compounds.add(v)
                    elif isinstance(value, str) and value.startswith('C'):
                        kegg_compounds.add(value)
        
        # Also check metabolite ID for KEGG patterns
        if metabolite.id.startswith('cpd'):
            # ModelSEED compound ID pattern
            continue
        elif re.match(r'C\d{5}', metabolite.id):
            kegg_compounds.add(metabolite.id)
    
    return kegg_compounds


def extract_kegg_reactions(model: cobra.Model) -> Set[str]:
    """Extract KEGG reaction IDs from reaction annotations."""
    kegg_reactions = set()
    
    for reaction in model.reactions:
        # Check annotations for KEGG reaction IDs
        if hasattr(reaction, 'annotation'):
            for key, value in reaction.annotation.items():
                if 'kegg.reaction' in key.lower():
                    if isinstance(value, list):
                        for v in value:
                            if isinstance(v, str) and v.startswith('R'):
                                kegg_reactions.add(v)
                    elif isinstance(value, str) and value.startswith('R'):
                        kegg_reactions.add(value)
        
        # Also check reaction ID for KEGG patterns
        if re.match(r'R\d{5}', reaction.id):
            kegg_reactions.add(reaction.id)
    
    return kegg_reactions


def extract_ec_numbers(model: cobra.Model) -> Set[str]:
    """Extract EC numbers from reaction annotations."""
    ec_numbers = set()
    
    for reaction in model.reactions:
        # Check annotations for EC numbers
        if hasattr(reaction, 'annotation'):
            for key, value in reaction.annotation.items():
                if 'ec-code' in key.lower() or 'enzyme' in key.lower():
                    if isinstance(value, list):
                        for v in value:
                            if isinstance(v, str) and re.match(r'\d+\.\d+\.\d+\.\d+', v):
                                ec_numbers.add(v)
                    elif isinstance(value, str) and re.match(r'\d+\.\d+\.\d+\.\d+', value):
                        ec_numbers.add(value)
        
        # Check reaction name/ID for EC patterns
        if hasattr(reaction, 'name') and reaction.name:
            ec_matches = re.findall(r'\d+\.\d+\.\d+\.\d+', reaction.name)
            ec_numbers.update(ec_matches)
    
    return ec_numbers


def extract_gene_names(model: cobra.Model) -> Set[str]:
    """Extract gene names for UniProt searches."""
    gene_names = set()
    
    for gene in model.genes:
        # Use gene ID as search term
        gene_names.add(gene.id)
        
        # Also check gene name if different from ID
        if hasattr(gene, 'name') and gene.name and gene.name != gene.id:
            gene_names.add(gene.name)
    
    return gene_names


def extract_compound_names(model: cobra.Model) -> Set[str]:
    """Extract metabolite names for ChEBI searches."""
    compound_names = set()
    
    for metabolite in model.metabolites:
        if hasattr(metabolite, 'name') and metabolite.name:
            # Clean up compound name
            name = metabolite.name.strip()
            if name and len(name) > 2:  # Skip very short names
                compound_names.add(name)
    
    return compound_names


def write_identifier_list(identifiers: Set[str], output_path: Path) -> None:
    """Write set of identifiers to file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with output_path.open('w') as f:
        f.write("# Generated identifier list\n")
        for identifier in sorted(identifiers):
            f.write(f"{identifier}\n")
    
    print(f"Wrote {len(identifiers)} identifiers to {output_path}")


def main() -> int:
    args = parse_args()
    
    print(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)
    
    output_dir = Path(args.output_dir)
    
    # Extract different types of identifiers
    print("Extracting identifiers...")
    
    kegg_compounds = extract_kegg_compounds(model)
    if kegg_compounds:
        write_identifier_list(kegg_compounds, output_dir / "kegg_compounds.txt")
    
    kegg_reactions = extract_kegg_reactions(model)
    if kegg_reactions:
        write_identifier_list(kegg_reactions, output_dir / "kegg_reactions.txt")
    
    ec_numbers = extract_ec_numbers(model)
    if ec_numbers:
        write_identifier_list(ec_numbers, output_dir / "ec_numbers.txt")
    
    gene_names = extract_gene_names(model)
    if gene_names:
        write_identifier_list(gene_names, output_dir / "gene_names.txt")
    
    compound_names = extract_compound_names(model)
    if compound_names:
        write_identifier_list(compound_names, output_dir / "compound_names.txt")
    
    print(f"\nSummary:")
    print(f"  KEGG compounds: {len(kegg_compounds)}")
    print(f"  KEGG reactions: {len(kegg_reactions)}")
    print(f"  EC numbers: {len(ec_numbers)}")
    print(f"  Gene names: {len(gene_names)}")
    print(f"  Compound names: {len(compound_names)}")
    
    print(f"\nExample API mapping commands:")
    if kegg_compounds:
        print(f"  python scripts/api_database_mapper.py kegg-compounds --input data/mappings/kegg_compounds.txt --output data/mappings/kegg_compounds_info.tsv")
    if kegg_reactions:
        print(f"  python scripts/api_database_mapper.py kegg-reactions --input data/mappings/kegg_reactions.txt --output data/mappings/kegg_reactions_info.tsv")
    if ec_numbers:
        print(f"  python scripts/api_database_mapper.py kegg-enzymes --input data/mappings/ec_numbers.txt --output data/mappings/ec_info.tsv")
    if gene_names:
        print(f"  python scripts/api_database_mapper.py uniprot-search --input data/mappings/gene_names.txt --output data/mappings/gene_uniprot_info.tsv --organism 'Solanum tuberosum'")
    if compound_names:
        print(f"  python scripts/api_database_mapper.py chebi-search --input data/mappings/compound_names.txt --output data/mappings/compound_chebi_info.tsv")
    
    return 0


if __name__ == "__main__":
    exit(main())