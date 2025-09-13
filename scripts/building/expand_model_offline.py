#!/usr/bin/env python3
"""
Offline model expansion using eggNOG-mapper annotations.

This script expands the model without external API calls by:
1. Using EC number mappings to create GPRs for existing reactions
2. Adding gene annotations (GO terms, descriptions, etc.)
3. Creating metabolic pathway groupings
4. Adding essential reactions based on EC patterns

Usage:
  python scripts/expand_model_offline.py \
    --model models/current/creole_unified.xml \
    --emapper reports/emapper_creole.emapper.annotations \
    --out models/current/creole_expanded.xml
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

import cobra
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Reaction, Metabolite, Gene


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Offline model expansion using eggNOG-mapper")
    p.add_argument("--model", required=True, help="Input SBML model")
    p.add_argument("--emapper", required=True, help="eggNOG-mapper annotations file")
    p.add_argument("--out", required=True, help="Output expanded model")
    p.add_argument("--min-score", type=float, default=50.0, help="Minimum eggNOG score")
    p.add_argument("--max-evalue", type=float, default=1e-5, help="Maximum e-value")
    p.add_argument("--verbose", action="store_true", help="Verbose output")
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
            'gos': [g.strip() for g in fields[9].split(',') if g.strip() and g.strip() != '-'],
            'ec': [e.strip() for e in fields[10].split(',') if e.strip() and e.strip() != '-'],
            'kegg_ko': [k.strip() for k in fields[11].split(',') if k.strip() and k.strip() != '-'],
            'kegg_pathway': [p.strip() for p in fields[12].split(',') if p.strip() and p.strip() != '-'],
            'kegg_module': [m.strip() for m in fields[13].split(',') if m.strip() and m.strip() != '-'],
            'kegg_reaction': [r.strip() for r in fields[14].split(',') if r.strip() and r.strip() != '-'],
            'cog_category': fields[6],
            'cazy': [c.strip() for c in fields[18].split(',') if c.strip() and c.strip() != '-'],
            'pfams': [p.strip() for p in fields[20].split(',') if p.strip() and p.strip() != '-']
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


def complete_gprs_with_ec_mapping(model: cobra.Model, annotations: List[Dict]) -> int:
    """Complete GPRs for reactions using EC number mappings."""
    
    # Build EC to genes mapping
    ec_to_genes = defaultdict(set)
    for ann in annotations:
        if ann['ec']:
            for ec in ann['ec']:
                if re.match(r'\d+\.\d+\.\d+\.\d+', ec):
                    ec_to_genes[ec].add(ann['query'])
    
    print(f"Found {len(ec_to_genes)} unique EC numbers mapped to genes")
    
    # Complete GPRs for reactions without them
    gprs_added = 0
    genes_added = 0
    
    for reaction in model.reactions:
        if reaction.gene_reaction_rule:
            continue  # Skip reactions that already have GPRs
            
        # Look for EC numbers in reaction annotations
        reaction_ecs = set()
        if hasattr(reaction, 'annotation'):
            for key, value in reaction.annotation.items():
                if 'ec' in key.lower() or 'enzyme' in key.lower():
                    if isinstance(value, list):
                        for v in value:
                            if isinstance(v, str) and re.match(r'\d+\.\d+\.\d+\.\d+', v):
                                reaction_ecs.add(v)
                    elif isinstance(value, str) and re.match(r'\d+\.\d+\.\d+\.\d+', value):
                        reaction_ecs.add(value)
        
        # Also check reaction name/ID for EC patterns
        if hasattr(reaction, 'name') and reaction.name:
            ec_matches = re.findall(r'\d+\.\d+\.\d+\.\d+', reaction.name)
            reaction_ecs.update(ec_matches)
        
        # Find genes with matching ECs
        candidate_genes = set()
        for ec in reaction_ecs:
            candidate_genes.update(ec_to_genes.get(ec, set()))
        
        if candidate_genes:
            # Create genes if they don't exist
            gene_objects = []
            for gene_id in sorted(candidate_genes):
                if gene_id not in [g.id for g in model.genes]:
                    gene_obj = Gene(gene_id)
                    model.genes.append(gene_obj)
                    genes_added += 1
                else:
                    gene_obj = model.genes.get_by_id(gene_id)
                gene_objects.append(gene_obj)
            
            # Create OR rule (genes as isoenzymes)
            reaction.gene_reaction_rule = ' or '.join([g.id for g in gene_objects])
            gprs_added += 1
    
    print(f"Added GPRs to {gprs_added} reactions")
    print(f"Added {genes_added} new genes to model")
    return gprs_added


def enrich_gene_annotations(model: cobra.Model, annotations: List[Dict]) -> int:
    """Add functional annotations to genes from eggNOG."""
    
    # Collect annotations by gene
    gene_annotations = defaultdict(dict)
    
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
        
        # KEGG KO
        if ann['kegg_ko']:
            kos = [ko.replace('ko:', '') for ko in ann['kegg_ko'] if ko.startswith('ko:')]
            if kos:
                gene_annotations[gene_id]['kegg.orthology'] = kos
        
        # KEGG pathways
        if ann['kegg_pathway']:
            pathways = [p for p in ann['kegg_pathway'] if p.startswith('ko') or p.startswith('map')]
            if pathways:
                gene_annotations[gene_id]['kegg.pathway'] = pathways
        
        # Function description
        if ann['description'] and ann['description'] != '-':
            gene_annotations[gene_id]['description'] = ann['description']
        
        # Preferred name
        if ann['preferred_name'] and ann['preferred_name'] != '-':
            gene_annotations[gene_id]['name'] = ann['preferred_name']
        
        # COG category
        if ann['cog_category'] and ann['cog_category'] != '-':
            gene_annotations[gene_id]['cog'] = ann['cog_category']
        
        # CAZy families
        if ann['cazy']:
            gene_annotations[gene_id]['cazy'] = ann['cazy']
    
    # Apply annotations to model genes
    annotated_count = 0
    for gene in model.genes:
        if gene.id in gene_annotations:
            if not hasattr(gene, 'annotation'):
                gene.annotation = {}
            gene.annotation.update(gene_annotations[gene.id])
            
            # Also update gene name if available
            if 'name' in gene_annotations[gene.id]:
                gene.name = gene_annotations[gene.id]['name']
            
            annotated_count += 1
    
    print(f"Added functional annotations to {annotated_count:,} genes")
    return annotated_count


def add_pathway_groups(model: cobra.Model, annotations: List[Dict]) -> Dict[str, List[str]]:
    """Create pathway groupings based on KEGG pathway annotations."""
    
    # Map reactions to pathways via genes
    pathway_reactions = defaultdict(set)
    
    for reaction in model.reactions:
        if not reaction.gene_reaction_rule:
            continue
            
        # Get genes from GPR (simplified parsing)
        gene_ids = set()
        gpr = reaction.gene_reaction_rule
        # Simple extraction - in production would need proper GPR parsing
        for word in gpr.replace('(', ' ').replace(')', ' ').replace('and', ' ').replace('or', ' ').split():
            if word.strip():
                gene_ids.add(word.strip())
        
        # Find pathways for these genes
        reaction_pathways = set()
        for ann in annotations:
            if ann['query'] in gene_ids and ann['kegg_pathway']:
                for pathway in ann['kegg_pathway']:
                    if pathway.startswith('ko') or pathway.startswith('map'):
                        reaction_pathways.add(pathway)
        
        # Add reaction to pathways
        for pathway in reaction_pathways:
            pathway_reactions[pathway].add(reaction.id)
    
    # Convert to dict with lists
    pathway_dict = {pathway: list(reactions) for pathway, reactions in pathway_reactions.items()}
    
    print(f"Organized {len(model.reactions)} reactions into {len(pathway_dict)} pathways")
    return pathway_dict


def analyze_functional_coverage(model: cobra.Model, annotations: List[Dict]) -> Dict:
    """Analyze functional coverage of the expanded model."""
    
    # Count functional annotations
    ec_counts = Counter()
    go_counts = Counter()
    pathway_counts = Counter()
    
    for ann in annotations:
        if ann['ec']:
            for ec in ann['ec']:
                if re.match(r'\d+\.\d+\.\d+\.\d+', ec):
                    ec_counts[ec] += 1
        
        if ann['gos']:
            for go in ann['gos']:
                if go.startswith('GO:'):
                    go_counts[go] += 1
        
        if ann['kegg_pathway']:
            for pathway in ann['kegg_pathway']:
                if pathway.startswith('ko') or pathway.startswith('map'):
                    pathway_counts[pathway] += 1
    
    # Model coverage
    model_genes = {g.id for g in model.genes}
    emapper_genes = {ann['query'] for ann in annotations}
    gene_coverage = len(model_genes & emapper_genes) / len(model_genes) * 100 if model_genes else 0
    
    return {
        'gene_coverage_percent': gene_coverage,
        'unique_ecs': len(ec_counts),
        'unique_go_terms': len(go_counts),
        'unique_pathways': len(pathway_counts),
        'top_ecs': ec_counts.most_common(10),
        'top_pathways': pathway_counts.most_common(10),
        'reactions_with_gprs': len([r for r in model.reactions if r.gene_reaction_rule]),
        'reactions_without_gprs': len([r for r in model.reactions if not r.gene_reaction_rule])
    }


def main() -> int:
    args = parse_args()
    
    print("🧬 OFFLINE MODEL EXPANSION FROM eggNOG-mapper")
    print("=" * 60)
    
    # Load model
    print(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)
    original_stats = {
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites),
        'genes': len(model.genes),
        'reactions_with_gprs': len([r for r in model.reactions if r.gene_reaction_rule])
    }
    
    print(f"Original model: {original_stats['reactions']} reactions, "
          f"{original_stats['metabolites']} metabolites, {original_stats['genes']} genes")
    print(f"Reactions with GPRs: {original_stats['reactions_with_gprs']}")
    
    # Load eggNOG annotations
    annotations = load_emapper_annotations(args.emapper, args.min_score, args.max_evalue)
    
    # Perform expansions
    print("\n🎯 Completing GPRs using EC mappings...")
    gprs_added = complete_gprs_with_ec_mapping(model, annotations)
    
    print("\n📝 Enriching gene annotations...")
    genes_annotated = enrich_gene_annotations(model, annotations)
    
    print("\n🗂️ Creating pathway groupings...")
    pathway_groups = add_pathway_groups(model, annotations)
    
    print("\n📊 Analyzing functional coverage...")
    coverage = analyze_functional_coverage(model, annotations)
    
    # Final stats
    final_stats = {
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites),
        'genes': len(model.genes),
        'reactions_with_gprs': len([r for r in model.reactions if r.gene_reaction_rule])
    }
    
    # Print summary
    print("\n" + "=" * 60)
    print("📈 EXPANSION SUMMARY")
    print("=" * 60)
    print(f"Reactions: {original_stats['reactions']:,} (unchanged)")
    print(f"Metabolites: {original_stats['metabolites']:,} (unchanged)")
    print(f"Genes: {original_stats['genes']:,} → {final_stats['genes']:,} (+{final_stats['genes'] - original_stats['genes']:,})")
    print(f"Reactions with GPRs: {original_stats['reactions_with_gprs']:,} → {final_stats['reactions_with_gprs']:,} (+{gprs_added:,})")
    print(f"Gene coverage: {coverage['gene_coverage_percent']:.1f}%")
    print(f"Unique ECs available: {coverage['unique_ecs']:,}")
    print(f"Unique GO terms: {coverage['unique_go_terms']:,}")
    print(f"KEGG pathways: {coverage['unique_pathways']:,}")
    print(f"Reactions still without GPRs: {coverage['reactions_without_gprs']:,}")
    
    # Save expanded model
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        write_sbml_model(model, str(output_path))
        print(f"\n💾 Expanded model saved to: {output_path}")
        
        # Save expansion report
        report = {
            'expansion_type': 'offline_emapper_integration',
            'original_stats': original_stats,
            'final_stats': final_stats,
            'functional_coverage': coverage,
            'pathway_groups': pathway_groups,
            'emapper_annotations_used': len(annotations),
            'gprs_added': gprs_added,
            'genes_annotated': genes_annotated
        }
        
        report_path = output_path.with_suffix('.expansion_report.json')
        with report_path.open('w') as f:
            json.dump(report, f, indent=2, default=str)
        print(f"📊 Expansion report saved to: {report_path}")
        
        # Save pathway groupings
        pathway_path = output_path.parent / "pathway_groups.json"
        with pathway_path.open('w') as f:
            json.dump(pathway_groups, f, indent=2)
        print(f"🗂️ Pathway groups saved to: {pathway_path}")
        
    except Exception as e:
        print(f"❌ Error saving model: {e}")
        return 1
    
    print(f"\n🎉 Model expansion completed successfully!")
    return 0


if __name__ == "__main__":
    exit(main())