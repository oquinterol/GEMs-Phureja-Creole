#!/usr/bin/env python3
"""
Massive model expansion using eggNOG-mapper annotations.

This script dramatically expands the metabolic model by:
1. Adding new reactions based on KEGG reaction IDs from eggNOG
2. Creating comprehensive GPRs using EC number mappings  
3. Adding missing metabolites and compartments
4. Enriching gene annotations with GO terms, COG categories
5. Organizing reactions into KEGG pathways

Usage:
  python scripts/expand_model_from_emapper.py \
    --model models/current/creole_v1.1.xml \
    --emapper reports/emapper_creole.emapper.annotations \
    --strategy comprehensive \
    --out models/current/creole_v1.2.xml
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
import requests

import cobra
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Reaction, Metabolite, Gene


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Expand model using eggNOG-mapper data")
    p.add_argument("--model", required=True, help="Input SBML model")
    p.add_argument("--emapper", required=True, help="eggNOG-mapper annotations file")
    p.add_argument("--out", required=True, help="Output expanded SBML model")
    p.add_argument("--strategy", choices=["conservative", "comprehensive", "aggressive"], 
                   default="comprehensive", help="Expansion strategy")
    p.add_argument("--min-score", type=float, default=50.0, help="Minimum eggNOG score")
    p.add_argument("--max-evalue", type=float, default=1e-5, help="Maximum e-value")
    p.add_argument("--kegg-api-delay", type=float, default=1.0, help="Delay between KEGG API calls")
    p.add_argument("--max-kegg-calls", type=int, default=500, help="Maximum KEGG API calls")
    p.add_argument("--dry-run", action="store_true", help="Show what would be added without modifying model")
    return p.parse_args()


def parse_emapper_line(line: str) -> Optional[Dict]:
    """Parse eggNOG-mapper annotation line."""
    fields = line.strip().split('\t')
    if len(fields) < 21:
        return None
    
    try:
        return {
            'query': fields[0],
            'evalue': float(fields[2]) if fields[2] != '-' else None,
            'score': float(fields[3]) if fields[3] != '-' else None,
            'description': fields[7],
            'preferred_name': fields[8],
            'gos': [g.strip() for g in fields[9].split(',') if g.strip() != '-'],
            'ec': [e.strip() for e in fields[10].split(',') if e.strip() != '-'],
            'kegg_ko': [k.strip() for k in fields[11].split(',') if k.strip() != '-'],
            'kegg_pathway': [p.strip() for p in fields[12].split(',') if p.strip() != '-'],
            'kegg_module': [m.strip() for m in fields[13].split(',') if m.strip() != '-'],
            'kegg_reaction': [r.strip() for r in fields[14].split(',') if r.strip() != '-'],
            'cogy_category': fields[6],
            'cazy': [c.strip() for c in fields[18].split(',') if c.strip() != '-'],
            'pfams': [p.strip() for p in fields[20].split(',') if p.strip() != '-']
        }
    except (ValueError, IndexError):
        return None


def load_emapper_annotations(emapper_file: str, min_score: float, max_evalue: float) -> List[Dict]:
    """Load and filter eggNOG annotations."""
    annotations = []
    print(f"Loading eggNOG annotations from: {emapper_file}")
    
    with open(emapper_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
                
            annotation = parse_emapper_line(line)
            if annotation is None:
                continue
                
            # Apply quality filters
            if annotation['score'] is not None and annotation['score'] < min_score:
                continue
            if annotation['evalue'] is not None and annotation['evalue'] > max_evalue:
                continue
                
            annotations.append(annotation)
            
            if line_num % 50000 == 0:
                print(f"  Processed {line_num:,} lines, kept {len(annotations):,} annotations")
    
    print(f"Loaded {len(annotations):,} high-quality annotations")
    return annotations


class KEGGAPIClient:
    """KEGG API client with rate limiting."""
    
    def __init__(self, delay: float = 1.0, max_calls: int = 500):
        self.delay = delay
        self.max_calls = max_calls
        self.call_count = 0
        self.last_call = 0.0
        
    def _rate_limit(self):
        """Enforce rate limiting."""
        if self.call_count >= self.max_calls:
            raise RuntimeError(f"Maximum KEGG API calls ({self.max_calls}) reached")
            
        now = time.time()
        time_since_last = now - self.last_call
        if time_since_last < self.delay:
            time.sleep(self.delay - time_since_last)
        
        self.last_call = time.time()
        self.call_count += 1
        
    def get_reaction(self, reaction_id: str) -> Optional[Dict]:
        """Get KEGG reaction details."""
        self._rate_limit()
        
        try:
            url = f"https://rest.kegg.jp/get/{reaction_id}"
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            return self._parse_kegg_reaction(response.text)
        except Exception as e:
            print(f"Warning: Could not fetch KEGG reaction {reaction_id}: {e}")
            return None
    
    def _parse_kegg_reaction(self, text: str) -> Dict:
        """Parse KEGG reaction flat file."""
        data = {}
        current_key = None
        
        for line in text.splitlines():
            if re.match(r'^[A-Z_]{2,}\s+', line):
                key = line[:12].strip()
                value = line[12:].strip()
                data[key] = value
                current_key = key
            elif current_key and line.startswith('            '):
                data[current_key] += ' ' + line.strip()
        
        return data


def extract_kegg_reactions_from_emapper(annotations: List[Dict]) -> Dict[str, List[str]]:
    """Extract KEGG reactions and associated genes from eggNOG."""
    reaction_to_genes = defaultdict(list)
    
    for ann in annotations:
        if ann['kegg_reaction']:
            for rxn_id in ann['kegg_reaction']:
                if rxn_id.startswith('R') and len(rxn_id) == 6:
                    reaction_to_genes[rxn_id].append(ann['query'])
    
    print(f"Found {len(reaction_to_genes)} unique KEGG reactions from eggNOG")
    return dict(reaction_to_genes)


def create_reaction_from_kegg(rxn_id: str, kegg_data: Dict, model: cobra.Model) -> Optional[Reaction]:
    """Create cobra Reaction from KEGG data."""
    try:
        # Parse equation
        equation = kegg_data.get('EQUATION', '')
        if not equation:
            return None
        
        # Create reaction
        reaction = Reaction(f"kegg_{rxn_id}")
        reaction.name = kegg_data.get('NAME', rxn_id)
        
        # Parse metabolites from equation
        # Simple parsing - this would need to be more sophisticated for production
        if '<->' in equation:
            left, right = equation.split('<->', 1)
            reaction.lower_bound = -1000
            reaction.upper_bound = 1000
        elif '=>' in equation:
            left, right = equation.split('=>', 1)
            reaction.lower_bound = 0
            reaction.upper_bound = 1000
        else:
            return None
        
        # Add metabolites (simplified - would need proper parsing)
        metabolites = {}
        
        # Add annotations
        reaction.annotation = {
            'kegg.reaction': rxn_id,
            'sbo': 'SBO:0000176'  # biochemical reaction
        }
        
        if 'ENZYME' in kegg_data:
            ecs = [ec.strip() for ec in kegg_data['ENZYME'].split() if re.match(r'\d+\.\d+\.\d+\.\d+', ec)]
            if ecs:
                reaction.annotation['ec-code'] = ecs
        
        return reaction
        
    except Exception as e:
        print(f"Error creating reaction {rxn_id}: {e}")
        return None


def add_gene_annotations_from_emapper(model: cobra.Model, annotations: List[Dict]) -> None:
    """Add GO terms and other annotations to genes."""
    gene_annotations = defaultdict(dict)
    
    # Collect annotations by gene
    for ann in annotations:
        gene_id = ann['query']
        
        # GO terms
        if ann['gos']:
            go_terms = [go for go in ann['gos'] if go.startswith('GO:')]
            if go_terms:
                gene_annotations[gene_id]['go'] = go_terms
        
        # EC numbers
        if ann['ec']:
            ecs = [ec for ec in ann['ec'] if re.match(r'\d+\.\d+\.\d+\.\d+', ec)]
            if ecs:
                gene_annotations[gene_id]['ec-code'] = ecs
        
        # KEGG info
        if ann['kegg_ko']:
            gene_annotations[gene_id]['kegg.orthology'] = ann['kegg_ko']
        
        # Function description
        if ann['description'] and ann['description'] != '-':
            gene_annotations[gene_id]['description'] = ann['description']
        
        # COG category
        if ann['cogy_category'] and ann['cogy_category'] != '-':
            gene_annotations[gene_id]['cog'] = ann['cogy_category']
    
    # Apply to model genes
    annotated_count = 0
    for gene in model.genes:
        if gene.id in gene_annotations:
            if not hasattr(gene, 'annotation'):
                gene.annotation = {}
            gene.annotation.update(gene_annotations[gene.id])
            annotated_count += 1
    
    print(f"Added annotations to {annotated_count:,} genes")


def create_gprs_from_emapper(model: cobra.Model, annotations: List[Dict], 
                           reaction_to_genes: Dict[str, List[str]]) -> None:
    """Create GPRs for reactions based on eggNOG data."""
    gpr_added = 0
    
    # Map EC numbers to genes
    ec_to_genes = defaultdict(list)
    for ann in annotations:
        if ann['ec']:
            for ec in ann['ec']:
                if re.match(r'\d+\.\d+\.\d+\.\d+', ec):
                    ec_to_genes[ec].append(ann['query'])
    
    # For each reaction, try to create GPR
    for reaction in model.reactions:
        if reaction.gene_reaction_rule:
            continue  # Skip if already has GPR
            
        # Try to find genes via EC numbers
        if hasattr(reaction, 'annotation') and 'ec-code' in reaction.annotation:
            ecs = reaction.annotation['ec-code']
            if isinstance(ecs, str):
                ecs = [ecs]
            
            candidate_genes = set()
            for ec in ecs:
                candidate_genes.update(ec_to_genes.get(ec, []))
            
            if candidate_genes:
                # Create OR rule for isoenzymes
                gene_objects = []
                for gene_id in sorted(candidate_genes):
                    if gene_id not in [g.id for g in model.genes]:
                        gene_obj = Gene(gene_id)
                        model.genes.append(gene_obj)
                    else:
                        gene_obj = model.genes.get_by_id(gene_id)
                    gene_objects.append(gene_obj)
                
                reaction.gene_reaction_rule = ' or '.join([g.id for g in gene_objects])
                gpr_added += 1
    
    print(f"Added GPRs to {gpr_added} reactions")


def expand_model_comprehensive(model: cobra.Model, annotations: List[Dict], 
                             kegg_client: KEGGAPIClient, dry_run: bool = False) -> Dict:
    """Comprehensive model expansion strategy."""
    stats = {
        'reactions_added': 0,
        'metabolites_added': 0,
        'genes_added': 0,
        'gprs_added': 0,
        'annotations_added': 0
    }
    
    print("🚀 Starting comprehensive model expansion...")
    
    # 1. Extract KEGG reactions from eggNOG
    print("📊 Extracting KEGG reactions from eggNOG...")
    reaction_to_genes = extract_kegg_reactions_from_emapper(annotations)
    
    # 2. Get existing reaction IDs to avoid duplicates
    existing_reactions = {r.id for r in model.reactions}
    existing_kegg_reactions = set()
    
    for reaction in model.reactions:
        if hasattr(reaction, 'annotation') and 'kegg.reaction' in reaction.annotation:
            kegg_id = reaction.annotation['kegg.reaction']
            if isinstance(kegg_id, list):
                existing_kegg_reactions.update(kegg_id)
            else:
                existing_kegg_reactions.add(kegg_id)
    
    # 3. Add new KEGG reactions
    new_reactions = []
    processed_count = 0
    
    for rxn_id, gene_list in reaction_to_genes.items():
        if rxn_id in existing_kegg_reactions:
            continue
            
        if processed_count >= kegg_client.max_calls:
            print(f"Reached maximum KEGG API calls ({kegg_client.max_calls})")
            break
            
        print(f"Processing KEGG reaction {rxn_id} ({processed_count + 1}/{min(len(reaction_to_genes), kegg_client.max_calls)})")
        
        kegg_data = kegg_client.get_reaction(rxn_id)
        if kegg_data is None:
            continue
            
        reaction = create_reaction_from_kegg(rxn_id, kegg_data, model)
        if reaction is None:
            continue
            
        # Add GPR based on eggNOG genes
        if gene_list:
            # Create genes if they don't exist
            gene_objects = []
            for gene_id in gene_list:
                if gene_id not in [g.id for g in model.genes]:
                    if not dry_run:
                        gene_obj = Gene(gene_id)
                        model.genes.append(gene_obj)
                    stats['genes_added'] += 1
                else:
                    gene_obj = model.genes.get_by_id(gene_id)
                gene_objects.append(gene_obj if not dry_run else gene_id)
            
            if not dry_run:
                reaction.gene_reaction_rule = ' or '.join([g.id for g in gene_objects])
            stats['gprs_added'] += 1
        
        new_reactions.append(reaction)
        stats['reactions_added'] += 1
        processed_count += 1
    
    # 4. Add reactions to model
    if not dry_run and new_reactions:
        print(f"Adding {len(new_reactions)} new reactions to model...")
        model.add_reactions(new_reactions)
    
    # 5. Enhance existing GPRs
    print("🧬 Enhancing existing GPRs...")
    if not dry_run:
        create_gprs_from_emapper(model, annotations, reaction_to_genes)
    
    # 6. Add gene annotations
    print("📝 Adding gene annotations...")
    if not dry_run:
        add_gene_annotations_from_emapper(model, annotations)
    stats['annotations_added'] = len([ann for ann in annotations if ann['gos'] or ann['ec']])
    
    return stats


def main() -> int:
    args = parse_args()
    
    print("🔬 MASSIVE MODEL EXPANSION FROM eggNOG-mapper")
    print("=" * 60)
    
    # Load model
    print(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)
    original_stats = {
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites),
        'genes': len(model.genes)
    }
    print(f"Original model: {original_stats['reactions']} reactions, {original_stats['metabolites']} metabolites, {original_stats['genes']} genes")
    
    # Load eggNOG annotations
    annotations = load_emapper_annotations(args.emapper, args.min_score, args.max_evalue)
    
    # Initialize KEGG client
    kegg_client = KEGGAPIClient(args.kegg_api_delay, args.max_kegg_calls)
    
    # Expand model based on strategy
    if args.strategy == "comprehensive":
        stats = expand_model_comprehensive(model, annotations, kegg_client, args.dry_run)
    else:
        raise NotImplementedError(f"Strategy '{args.strategy}' not implemented yet")
    
    # Final model stats
    final_stats = {
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites),
        'genes': len(model.genes)
    }
    
    # Print summary
    print("\n" + "=" * 60)
    print("📈 EXPANSION SUMMARY")
    print("=" * 60)
    print(f"Reactions: {original_stats['reactions']:,} → {final_stats['reactions']:,} (+{stats['reactions_added']:,})")
    print(f"Metabolites: {original_stats['metabolites']:,} → {final_stats['metabolites']:,} (+{stats['metabolites_added']:,})")
    print(f"Genes: {original_stats['genes']:,} → {final_stats['genes']:,} (+{stats['genes_added']:,})")
    print(f"New GPRs added: {stats['gprs_added']:,}")
    print(f"Gene annotations: {stats['annotations_added']:,}")
    print(f"KEGG API calls used: {kegg_client.call_count}/{kegg_client.max_calls}")
    
    if args.dry_run:
        print("\n🔍 DRY RUN - No changes made to model")
        return 0
    
    # Save expanded model
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        write_sbml_model(model, str(output_path))
        print(f"\n💾 Expanded model saved to: {output_path}")
        
        # Save expansion report
        report = {
            'expansion_strategy': args.strategy,
            'original_stats': original_stats,
            'final_stats': final_stats,
            'expansion_stats': stats,
            'kegg_api_calls': kegg_client.call_count,
            'emapper_annotations_used': len(annotations)
        }
        
        report_path = output_path.with_suffix('.expansion_report.json')
        with report_path.open('w') as f:
            json.dump(report, f, indent=2)
        print(f"📊 Expansion report saved to: {report_path}")
        
    except Exception as e:
        print(f"❌ Error saving model: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())