#!/usr/bin/env python3
"""
API Database Mapper - Fetch and map metabolic information from various databases.

Supports multiple database APIs:
- KEGG: Compounds, reactions, pathways, enzymes
- ChEBI: Chemical entities, structures, ontology
- UniProt: Protein information, GO terms, EC numbers
- Rhea: Reaction information and cross-references
- BiGG: Metabolic model components

Usage:
  # Map KEGG compounds to ChEBI IDs
  python scripts/api_database_mapper.py kegg-compounds --input data/mappings/kegg_compounds.txt --output data/mappings/kegg2chebi.tsv
  
  # Get UniProt info for gene list
  python scripts/api_database_mapper.py uniprot-genes --input data/mappings/gene_list.txt --output data/mappings/gene_info.tsv
  
  # Map reactions to pathways
  python scripts/api_database_mapper.py kegg-pathways --input data/mappings/rxn_list.txt --output data/mappings/rxn2pathway.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Set, Union
import requests
from urllib.parse import urljoin
import xml.etree.ElementTree as ET


class RateLimiter:
    """Simple rate limiter to be respectful to APIs."""
    
    def __init__(self, calls_per_second: float = 1.0):
        self.calls_per_second = calls_per_second
        self.last_call = 0.0
    
    def wait(self):
        """Wait if necessary to respect rate limit."""
        now = time.time()
        time_since_last = now - self.last_call
        min_interval = 1.0 / self.calls_per_second
        
        if time_since_last < min_interval:
            time.sleep(min_interval - time_since_last)
        
        self.last_call = time.time()


class KeggAPI:
    """KEGG REST API client."""
    
    def __init__(self, rate_limit: float = 1.0):
        self.base_url = "http://rest.kegg.jp/"
        self.rate_limiter = RateLimiter(rate_limit)
    
    def _get(self, endpoint: str) -> Optional[str]:
        """Make GET request with rate limiting."""
        self.rate_limiter.wait()
        
        try:
            response = requests.get(urljoin(self.base_url, endpoint), timeout=30)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            print(f"KEGG API error for {endpoint}: {e}")
            return None
    
    def get_compound(self, compound_id: str) -> Optional[Dict]:
        """Get compound information."""
        data = self._get(f"get/cpd:{compound_id}")
        if not data:
            return None
        
        # Parse KEGG flat file format
        info = {"id": compound_id}
        lines = data.strip().split('\n')
        current_field = None
        
        for line in lines:
            if line.startswith('NAME'):
                info['name'] = line[12:].strip()
            elif line.startswith('FORMULA'):
                info['formula'] = line[12:].strip()
            elif line.startswith('DBLINKS'):
                current_field = 'dblinks'
                info['dblinks'] = {}
                dblink = line[12:].strip()
                if ':' in dblink:
                    db, ids = dblink.split(':', 1)
                    info['dblinks'][db.strip()] = [id.strip() for id in ids.split()]
            elif current_field == 'dblinks' and line.startswith('           '):
                dblink = line.strip()
                if ':' in dblink:
                    db, ids = dblink.split(':', 1)
                    info['dblinks'][db.strip()] = [id.strip() for id in ids.split()]
        
        return info
    
    def get_reaction(self, reaction_id: str) -> Optional[Dict]:
        """Get reaction information."""
        data = self._get(f"get/rn:{reaction_id}")
        if not data:
            return None
        
        info = {"id": reaction_id}
        lines = data.strip().split('\n')
        current_field = None
        
        for line in lines:
            if line.startswith('NAME'):
                info['name'] = line[12:].strip()
            elif line.startswith('DEFINITION'):
                info['definition'] = line[12:].strip()
            elif line.startswith('ENZYME'):
                info['enzyme'] = line[12:].strip().split()
            elif line.startswith('PATHWAY'):
                current_field = 'pathway'
                info['pathways'] = [line[12:].strip()]
            elif current_field == 'pathway' and line.startswith('           '):
                info['pathways'].append(line.strip())
            elif line.startswith('DBLINKS'):
                current_field = 'dblinks'
                info['dblinks'] = {}
                dblink = line[12:].strip()
                if ':' in dblink:
                    db, ids = dblink.split(':', 1)
                    info['dblinks'][db.strip()] = [id.strip() for id in ids.split()]
        
        return info
    
    def get_enzyme(self, ec_number: str) -> Optional[Dict]:
        """Get enzyme information."""
        data = self._get(f"get/ec:{ec_number}")
        if not data:
            return None
        
        info = {"ec": ec_number}
        lines = data.strip().split('\n')
        
        for line in lines:
            if line.startswith('NAME'):
                info['name'] = line[12:].strip()
            elif line.startswith('REACTION'):
                info['reactions'] = line[12:].strip().split()
        
        return info


class ChEBIAPI:
    """ChEBI API client."""
    
    def __init__(self, rate_limit: float = 1.0):
        self.base_url = "https://www.ebi.ac.uk/webservices/chebi/2.0/test/"
        self.rate_limiter = RateLimiter(rate_limit)
    
    def _get(self, endpoint: str, params: Dict = None) -> Optional[requests.Response]:
        """Make GET request with rate limiting."""
        self.rate_limiter.wait()
        
        try:
            response = requests.get(urljoin(self.base_url, endpoint), params=params, timeout=30)
            response.raise_for_status()
            return response
        except requests.RequestException as e:
            print(f"ChEBI API error for {endpoint}: {e}")
            return None
    
    def search_compound(self, query: str, search_category: str = "ALL") -> Optional[List[Dict]]:
        """Search for compounds."""
        params = {
            "query": query,
            "searchCategory": search_category,
            "maximumResults": 10
        }
        
        response = self._get("getLiteEntity", params)
        if not response:
            return None
        
        # Parse XML response
        try:
            root = ET.fromstring(response.text)
            entities = []
            
            for entity in root.findall(".//return"):
                chebi_id = entity.find("chebiId")
                chebi_name = entity.find("chebiAsciiName")
                
                if chebi_id is not None and chebi_name is not None:
                    entities.append({
                        "chebi_id": chebi_id.text,
                        "name": chebi_name.text
                    })
            
            return entities
            
        except ET.ParseError as e:
            print(f"ChEBI XML parse error: {e}")
            return None


class UniProtAPI:
    """UniProt API client."""
    
    def __init__(self, rate_limit: float = 1.0):
        self.base_url = "https://rest.uniprot.org/"
        self.rate_limiter = RateLimiter(rate_limit)
    
    def _get(self, endpoint: str, params: Dict = None) -> Optional[requests.Response]:
        """Make GET request with rate limiting."""
        self.rate_limiter.wait()
        
        try:
            response = requests.get(urljoin(self.base_url, endpoint), params=params, timeout=30)
            response.raise_for_status()
            return response
        except requests.RequestException as e:
            print(f"UniProt API error for {endpoint}: {e}")
            return None
    
    def search_proteins(self, query: str, organism: str = None) -> Optional[List[Dict]]:
        """Search for proteins."""
        search_query = query
        if organism:
            search_query += f" AND organism_name:{organism}"
        
        params = {
            "query": search_query,
            "format": "json",
            "size": 10
        }
        
        response = self._get("uniprotkb/search", params)
        if not response:
            return None
        
        try:
            data = response.json()
            return data.get("results", [])
        except json.JSONDecodeError as e:
            print(f"UniProt JSON parse error: {e}")
            return None


class RheaAPI:
    """Rhea reaction database API client."""
    
    def __init__(self, rate_limit: float = 1.0):
        self.base_url = "https://www.rhea-db.org/rhea/"
        self.rate_limiter = RateLimiter(rate_limit)
    
    def _get(self, endpoint: str, params: Dict = None) -> Optional[requests.Response]:
        """Make GET request with rate limiting."""
        self.rate_limiter.wait()
        
        try:
            response = requests.get(urljoin(self.base_url, endpoint), params=params, timeout=30)
            response.raise_for_status()
            return response
        except requests.RequestException as e:
            print(f"Rhea API error for {endpoint}: {e}")
            return None
    
    def get_reaction(self, rhea_id: str) -> Optional[Dict]:
        """Get reaction information."""
        response = self._get(f"{rhea_id}?format=json")
        if not response:
            return None
        
        try:
            return response.json()
        except json.JSONDecodeError as e:
            print(f"Rhea JSON parse error: {e}")
            return None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fetch and map information from metabolic databases")
    
    subparsers = p.add_subparsers(dest="command", help="Database mapping commands")
    
    # KEGG commands
    kegg_compounds = subparsers.add_parser("kegg-compounds", help="Map KEGG compounds to other databases")
    kegg_compounds.add_argument("--input", required=True, help="Input file with KEGG compound IDs")
    kegg_compounds.add_argument("--output", required=True, help="Output TSV file")
    kegg_compounds.add_argument("--rate-limit", type=float, default=1.0, help="Requests per second")
    
    kegg_reactions = subparsers.add_parser("kegg-reactions", help="Map KEGG reactions")
    kegg_reactions.add_argument("--input", required=True, help="Input file with KEGG reaction IDs")
    kegg_reactions.add_argument("--output", required=True, help="Output TSV file")
    kegg_reactions.add_argument("--rate-limit", type=float, default=1.0, help="Requests per second")
    
    kegg_enzymes = subparsers.add_parser("kegg-enzymes", help="Map EC numbers via KEGG")
    kegg_enzymes.add_argument("--input", required=True, help="Input file with EC numbers")
    kegg_enzymes.add_argument("--output", required=True, help="Output TSV file")
    kegg_enzymes.add_argument("--rate-limit", type=float, default=1.0, help="Requests per second")
    
    # ChEBI commands
    chebi_search = subparsers.add_parser("chebi-search", help="Search ChEBI for compound names")
    chebi_search.add_argument("--input", required=True, help="Input file with compound names")
    chebi_search.add_argument("--output", required=True, help="Output TSV file")
    chebi_search.add_argument("--rate-limit", type=float, default=1.0, help="Requests per second")
    
    # UniProt commands
    uniprot_search = subparsers.add_parser("uniprot-search", help="Search UniProt for genes/proteins")
    uniprot_search.add_argument("--input", required=True, help="Input file with gene names")
    uniprot_search.add_argument("--output", required=True, help="Output TSV file")
    uniprot_search.add_argument("--organism", help="Organism name filter")
    uniprot_search.add_argument("--rate-limit", type=float, default=1.0, help="Requests per second")
    
    return p.parse_args()


def read_input_file(file_path: str) -> List[str]:
    """Read input file and return list of identifiers."""
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    
    identifiers = []
    with path.open('r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                identifiers.append(line)
    
    return identifiers


def cmd_kegg_compounds(args) -> None:
    """Map KEGG compounds to other databases."""
    kegg = KeggAPI(args.rate_limit)
    compound_ids = read_input_file(args.input)
    
    results = []
    for compound_id in compound_ids:
        print(f"Processing compound: {compound_id}")
        info = kegg.get_compound(compound_id)
        
        if info:
            row = {
                "kegg_id": compound_id,
                "name": info.get("name", ""),
                "formula": info.get("formula", ""),
                "chebi_ids": "",
                "pubchem_ids": "",
                "cas_ids": ""
            }
            
            if "dblinks" in info:
                dblinks = info["dblinks"]
                if "ChEBI" in dblinks:
                    row["chebi_ids"] = ";".join(dblinks["ChEBI"])
                if "PubChem" in dblinks:
                    row["pubchem_ids"] = ";".join(dblinks["PubChem"])
                if "CAS" in dblinks:
                    row["cas_ids"] = ";".join(dblinks["CAS"])
            
            results.append(row)
        else:
            print(f"  Could not retrieve info for {compound_id}")
    
    # Write results
    with open(args.output, 'w', newline='') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
    
    print(f"Results written to: {args.output}")


def cmd_kegg_reactions(args) -> None:
    """Map KEGG reactions."""
    kegg = KeggAPI(args.rate_limit)
    reaction_ids = read_input_file(args.input)
    
    results = []
    for reaction_id in reaction_ids:
        print(f"Processing reaction: {reaction_id}")
        info = kegg.get_reaction(reaction_id)
        
        if info:
            row = {
                "kegg_id": reaction_id,
                "name": info.get("name", ""),
                "definition": info.get("definition", ""),
                "enzymes": ";".join(info.get("enzyme", [])),
                "pathways": ";".join(info.get("pathways", [])),
                "rhea_ids": ""
            }
            
            if "dblinks" in info and "Rhea" in info["dblinks"]:
                row["rhea_ids"] = ";".join(info["dblinks"]["Rhea"])
            
            results.append(row)
        else:
            print(f"  Could not retrieve info for {reaction_id}")
    
    # Write results
    with open(args.output, 'w', newline='') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
    
    print(f"Results written to: {args.output}")


def cmd_kegg_enzymes(args) -> None:
    """Map EC numbers via KEGG."""
    kegg = KeggAPI(args.rate_limit)
    ec_numbers = read_input_file(args.input)
    
    results = []
    for ec_number in ec_numbers:
        print(f"Processing EC: {ec_number}")
        info = kegg.get_enzyme(ec_number)
        
        if info:
            row = {
                "ec_number": ec_number,
                "name": info.get("name", ""),
                "reactions": ";".join(info.get("reactions", []))
            }
            results.append(row)
        else:
            print(f"  Could not retrieve info for {ec_number}")
    
    # Write results
    with open(args.output, 'w', newline='') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
    
    print(f"Results written to: {args.output}")


def cmd_chebi_search(args) -> None:
    """Search ChEBI for compound names."""
    chebi = ChEBIAPI(args.rate_limit)
    compound_names = read_input_file(args.input)
    
    results = []
    for compound_name in compound_names:
        print(f"Searching ChEBI for: {compound_name}")
        entities = chebi.search_compound(compound_name)
        
        if entities:
            for entity in entities:
                row = {
                    "query": compound_name,
                    "chebi_id": entity["chebi_id"],
                    "name": entity["name"]
                }
                results.append(row)
        else:
            print(f"  No results for {compound_name}")
    
    # Write results
    with open(args.output, 'w', newline='') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
    
    print(f"Results written to: {args.output}")


def cmd_uniprot_search(args) -> None:
    """Search UniProt for genes/proteins."""
    uniprot = UniProtAPI(args.rate_limit)
    gene_names = read_input_file(args.input)
    
    results = []
    for gene_name in gene_names:
        print(f"Searching UniProt for: {gene_name}")
        proteins = uniprot.search_proteins(gene_name, args.organism)
        
        if proteins:
            for protein in proteins[:3]:  # Limit to top 3 results
                row = {
                    "query": gene_name,
                    "uniprot_id": protein.get("primaryAccession", ""),
                    "gene_names": ";".join([g.get("geneName", {}).get("value", "") for g in protein.get("genes", [])]),
                    "protein_name": protein.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
                    "organism": protein.get("organism", {}).get("scientificName", ""),
                    "ec_numbers": ""
                }
                
                # Extract EC numbers
                ec_numbers = []
                for comment in protein.get("comments", []):
                    if comment.get("commentType") == "CATALYTIC_ACTIVITY":
                        ec = comment.get("reaction", {}).get("ecNumber")
                        if ec:
                            ec_numbers.append(ec)
                
                row["ec_numbers"] = ";".join(ec_numbers)
                results.append(row)
        else:
            print(f"  No results for {gene_name}")
    
    # Write results
    with open(args.output, 'w', newline='') as f:
        if results:
            writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
    
    print(f"Results written to: {args.output}")


def main() -> int:
    args = parse_args()
    
    if not args.command:
        print("Please specify a command. Use --help for options.")
        return 1
    
    try:
        if args.command == "kegg-compounds":
            cmd_kegg_compounds(args)
        elif args.command == "kegg-reactions":
            cmd_kegg_reactions(args)
        elif args.command == "kegg-enzymes":
            cmd_kegg_enzymes(args)
        elif args.command == "chebi-search":
            cmd_chebi_search(args)
        elif args.command == "uniprot-search":
            cmd_uniprot_search(args)
        else:
            print(f"Unknown command: {args.command}")
            return 1
            
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":
    exit(main())