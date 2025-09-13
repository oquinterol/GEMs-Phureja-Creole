#!/usr/bin/env python3
"""
Analyze eggNOG-mapper coverage and potential for model expansion.

This script performs comprehensive analysis of eggNOG-mapper results to understand:
- How many new reactions could be added
- Coverage of existing model reactions
- Potential for completing missing GPRs
- Functional annotation enrichment opportunities

Usage:
  python scripts/analyze_emapper_coverage.py \
    --emapper reports/emapper_creole.emapper.annotations \
    --model models/current/creole_unified.xml \
    --out reports/emapper_analysis.json
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple
import cobra
from cobra.io import read_sbml_model


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Analyze eggNOG-mapper coverage for model expansion")
    p.add_argument("--emapper", required=True, help="eggNOG-mapper annotations file (.emapper.annotations)")
    p.add_argument("--model", required=True, help="Current SBML model")
    p.add_argument("--out", required=True, help="Output analysis JSON file")
    p.add_argument("--min-score", type=float, default=50.0, help="Minimum eggNOG score threshold")
    p.add_argument("--max-evalue", type=float, default=1e-5, help="Maximum e-value threshold")
    return p.parse_args()


def parse_emapper_line(line: str) -> Dict:
    """Parse a single eggNOG-mapper annotation line."""
    fields = line.strip().split('\t')
    if len(fields) < 21:
        return None
    
    return {
        'query': fields[0],
        'seed_ortholog': fields[1],
        'evalue': float(fields[2]) if fields[2] != '-' else None,
        'score': float(fields[3]) if fields[3] != '-' else None,
        'eggnog_ogs': fields[4],
        'max_annot_lvl': fields[5],
        'cog_category': fields[6],
        'description': fields[7],
        'preferred_name': fields[8],
        'gos': fields[9].split(',') if fields[9] != '-' else [],
        'ec': fields[10].split(',') if fields[10] != '-' else [],
        'kegg_ko': fields[11].split(',') if fields[11] != '-' else [],
        'kegg_pathway': fields[12].split(',') if fields[12] != '-' else [],
        'kegg_module': fields[13].split(',') if fields[13] != '-' else [],
        'kegg_reaction': fields[14].split(',') if fields[14] != '-' else [],
        'kegg_rclass': fields[15].split(',') if fields[15] != '-' else [],
        'brite': fields[16].split(',') if fields[16] != '-' else [],
        'kegg_tc': fields[17].split(',') if fields[17] != '-' else [],
        'cazy': fields[18].split(',') if fields[18] != '-' else [],
        'bigg_reaction': fields[19].split(',') if fields[19] != '-' else [],
        'pfams': fields[20].split(',') if fields[20] != '-' else []
    }


def load_emapper_annotations(emapper_file: str, min_score: float = 50.0, max_evalue: float = 1e-5) -> List[Dict]:
    """Load and filter eggNOG-mapper annotations."""
    annotations = []
    total_lines = 0
    filtered_lines = 0
    
    print(f"Loading eggNOG annotations from: {emapper_file}")
    
    with open(emapper_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            total_lines += 1
            annotation = parse_emapper_line(line)
            
            if annotation is None:
                continue
                
            # Apply quality filters
            if annotation['score'] is not None and annotation['score'] < min_score:
                continue
            if annotation['evalue'] is not None and annotation['evalue'] > max_evalue:
                continue
            
            filtered_lines += 1
            annotations.append(annotation)
    
    print(f"Loaded {filtered_lines:,} high-quality annotations from {total_lines:,} total lines")
    return annotations


def analyze_model_coverage(model: cobra.Model, annotations: List[Dict]) -> Dict:
    """Analyze how well current model is covered by eggNOG annotations."""
    
    # Get model statistics
    model_stats = {
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites), 
        'genes': len(model.genes)
    }
    
    # Get model gene IDs
    model_gene_ids = {g.id for g in model.genes}
    
    # Get eggNOG gene IDs
    emapper_gene_ids = {ann['query'] for ann in annotations}
    
    # Calculate gene coverage
    genes_in_common = model_gene_ids & emapper_gene_ids
    genes_only_model = model_gene_ids - emapper_gene_ids
    genes_only_emapper = emapper_gene_ids - model_gene_ids
    
    gene_coverage = {
        'model_genes': len(model_gene_ids),
        'emapper_genes': len(emapper_gene_ids),
        'genes_in_common': len(genes_in_common),
        'genes_only_in_model': len(genes_only_model),
        'genes_only_in_emapper': len(genes_only_emapper),
        'coverage_percentage': len(genes_in_common) / len(model_gene_ids) * 100 if model_gene_ids else 0
    }
    
    # Analyze reactions without GPRs
    reactions_without_gprs = [r for r in model.reactions if not r.gene_reaction_rule]
    
    return {
        'model_stats': model_stats,
        'gene_coverage': gene_coverage,
        'reactions_without_gprs': len(reactions_without_gprs),
        'reactions_without_gprs_list': [r.id for r in reactions_without_gprs]
    }


def analyze_functional_content(annotations: List[Dict]) -> Dict:
    """Analyze functional content available in eggNOG annotations."""
    
    # Count functional annotations
    ec_counter = Counter()
    kegg_ko_counter = Counter()
    kegg_pathway_counter = Counter()
    kegg_reaction_counter = Counter()
    go_counter = Counter()
    cog_counter = Counter()
    cazy_counter = Counter()
    
    genes_with_ec = 0
    genes_with_kegg_ko = 0
    genes_with_kegg_reaction = 0
    genes_with_go = 0
    
    for ann in annotations:
        # EC numbers
        if ann['ec'] and ann['ec'] != ['-']:
            genes_with_ec += 1
            for ec in ann['ec']:
                if ec != '-' and re.match(r'\d+\.\d+\.\d+\.\d+', ec):
                    ec_counter[ec] += 1
        
        # KEGG KOs
        if ann['kegg_ko'] and ann['kegg_ko'] != ['-']:
            genes_with_kegg_ko += 1
            for ko in ann['kegg_ko']:
                if ko.startswith('ko:'):
                    kegg_ko_counter[ko] += 1
        
        # KEGG pathways
        for pathway in ann['kegg_pathway']:
            if pathway.startswith('ko') or pathway.startswith('map'):
                kegg_pathway_counter[pathway] += 1
        
        # KEGG reactions
        if ann['kegg_reaction'] and ann['kegg_reaction'] != ['-']:
            genes_with_kegg_reaction += 1
            for rxn in ann['kegg_reaction']:
                if rxn.startswith('R'):
                    kegg_reaction_counter[rxn] += 1
        
        # GO terms
        if ann['gos'] and ann['gos'] != ['-']:
            genes_with_go += 1
            for go_term in ann['gos']:
                if go_term.startswith('GO:'):
                    go_counter[go_term] += 1
        
        # COG categories
        if ann['cog_category'] and ann['cog_category'] != '-':
            cog_counter[ann['cog_category']] += 1
        
        # CAZy families
        for cazy in ann['cazy']:
            if cazy != '-':
                cazy_counter[cazy] += 1
    
    return {
        'total_genes': len(annotations),
        'ec_numbers': {
            'unique_ecs': len(ec_counter),
            'genes_with_ec': genes_with_ec,
            'total_ec_assignments': sum(ec_counter.values()),
            'top_10_ecs': ec_counter.most_common(10)
        },
        'kegg_ko': {
            'unique_kos': len(kegg_ko_counter),
            'genes_with_ko': genes_with_kegg_ko,
            'total_ko_assignments': sum(kegg_ko_counter.values())
        },
        'kegg_pathways': {
            'unique_pathways': len(kegg_pathway_counter),
            'total_pathway_assignments': sum(kegg_pathway_counter.values()),
            'top_10_pathways': kegg_pathway_counter.most_common(10)
        },
        'kegg_reactions': {
            'unique_reactions': len(kegg_reaction_counter),
            'genes_with_kegg_reaction': genes_with_kegg_reaction,
            'total_reaction_assignments': sum(kegg_reaction_counter.values()),
            'top_10_reactions': kegg_reaction_counter.most_common(10)
        },
        'go_terms': {
            'unique_go_terms': len(go_counter),
            'genes_with_go': genes_with_go,
            'total_go_assignments': sum(go_counter.values())
        },
        'cog_categories': {
            'unique_cog_categories': len(cog_counter),
            'total_cog_assignments': sum(cog_counter.values()),
            'cog_distribution': dict(cog_counter.most_common())
        },
        'cazy_families': {
            'unique_cazy_families': len(cazy_counter),
            'total_cazy_assignments': sum(cazy_counter.values())
        }
    }


def analyze_expansion_potential(model: cobra.Model, annotations: List[Dict]) -> Dict:
    """Analyze potential for model expansion using eggNOG data."""
    
    # Get current model ECs
    current_ecs = set()
    for reaction in model.reactions:
        if hasattr(reaction, 'annotation'):
            for key, value in reaction.annotation.items():
                if 'ec-code' in key.lower():
                    if isinstance(value, list):
                        current_ecs.update(value)
                    else:
                        current_ecs.add(value)
    
    # Get eggNOG ECs
    emapper_ecs = set()
    for ann in annotations:
        for ec in ann['ec']:
            if ec != '-' and re.match(r'\d+\.\d+\.\d+\.\d+', ec):
                emapper_ecs.add(ec)
    
    # Analyze EC coverage
    ec_overlap = current_ecs & emapper_ecs
    new_ecs = emapper_ecs - current_ecs
    missing_ecs = current_ecs - emapper_ecs
    
    # Get KEGG reactions from eggNOG
    kegg_reactions = set()
    for ann in annotations:
        for rxn in ann['kegg_reaction']:
            if rxn.startswith('R') and len(rxn) == 6:
                kegg_reactions.add(rxn)
    
    # Analyze potential GPR completion
    reactions_without_gprs = [r for r in model.reactions if not r.gene_reaction_rule]
    
    # For each reaction without GPR, check if we have genes with matching EC
    potential_gpr_completions = []
    for reaction in reactions_without_gprs:
        reaction_ecs = set()
        if hasattr(reaction, 'annotation'):
            for key, value in reaction.annotation.items():
                if 'ec-code' in key.lower():
                    if isinstance(value, list):
                        reaction_ecs.update(value)
                    else:
                        reaction_ecs.add(value)
        
        # Find genes with matching ECs
        matching_genes = []
        for ann in annotations:
            if any(ec in reaction_ecs for ec in ann['ec'] if ec != '-'):
                matching_genes.append(ann['query'])
        
        if matching_genes:
            potential_gpr_completions.append({
                'reaction_id': reaction.id,
                'reaction_ecs': list(reaction_ecs),
                'matching_genes': matching_genes[:10]  # Limit to first 10
            })
    
    return {
        'ec_analysis': {
            'current_model_ecs': len(current_ecs),
            'emapper_ecs': len(emapper_ecs),
            'ec_overlap': len(ec_overlap),
            'new_ecs_from_emapper': len(new_ecs),
            'model_ecs_not_in_emapper': len(missing_ecs),
            'coverage_percentage': len(ec_overlap) / len(current_ecs) * 100 if current_ecs else 0
        },
        'kegg_reactions': {
            'unique_kegg_reactions': len(kegg_reactions),
            'kegg_reactions_list': sorted(list(kegg_reactions))
        },
        'gpr_completion_potential': {
            'reactions_without_gprs': len(reactions_without_gprs),
            'reactions_with_potential_gprs': len(potential_gpr_completions),
            'completion_percentage': len(potential_gpr_completions) / len(reactions_without_gprs) * 100 if reactions_without_gprs else 0,
            'examples': potential_gpr_completions[:20]  # Show first 20 examples
        }
    }


def generate_expansion_recommendations(analysis: Dict) -> Dict:
    """Generate recommendations for model expansion based on analysis."""
    
    recommendations = {
        'high_priority': [],
        'medium_priority': [],
        'low_priority': []
    }
    
    # High priority recommendations
    gpr_completion = analysis['expansion_potential']['gpr_completion_potential']
    if gpr_completion['completion_percentage'] > 50:
        recommendations['high_priority'].append({
            'action': 'Complete GPRs using eggNOG EC mappings',
            'impact': f"Could complete GPRs for {gpr_completion['reactions_with_potential_gprs']} reactions ({gpr_completion['completion_percentage']:.1f}%)",
            'difficulty': 'Medium'
        })
    
    # Check for significant functional content
    func_content = analysis['functional_content']
    if func_content['ec_numbers']['unique_ecs'] > 1000:
        recommendations['high_priority'].append({
            'action': 'Add new reactions based on eggNOG EC numbers',
            'impact': f"Potential to add reactions from {func_content['ec_numbers']['unique_ecs']} unique ECs",
            'difficulty': 'High'
        })
    
    if func_content['kegg_reactions']['unique_reactions'] > 100:
        recommendations['high_priority'].append({
            'action': 'Add KEGG reactions from eggNOG annotations',
            'impact': f"Potential to add {func_content['kegg_reactions']['unique_reactions']} KEGG reactions",
            'difficulty': 'High'
        })
    
    # Medium priority
    if func_content['go_terms']['genes_with_go'] > 1000:
        recommendations['medium_priority'].append({
            'action': 'Enrich gene annotations with GO terms',
            'impact': f"Add GO terms to {func_content['go_terms']['genes_with_go']} genes",
            'difficulty': 'Low'
        })
    
    if func_content['kegg_pathways']['unique_pathways'] > 50:
        recommendations['medium_priority'].append({
            'action': 'Add pathway information to model',
            'impact': f"Organize reactions into {func_content['kegg_pathways']['unique_pathways']} pathways",
            'difficulty': 'Medium'
        })
    
    return recommendations


def main() -> int:
    args = parse_args()
    
    print("🔍 Analyzing eggNOG-mapper coverage for model expansion...")
    print("=" * 60)
    
    # Load current model
    print(f"Loading current model: {args.model}")
    model = read_sbml_model(args.model)
    print(f"Model stats: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, {len(model.genes)} genes")
    
    # Load eggNOG annotations
    annotations = load_emapper_annotations(args.emapper, args.min_score, args.max_evalue)
    
    # Perform analyses
    print("\n📊 Analyzing model coverage...")
    model_coverage = analyze_model_coverage(model, annotations)
    
    print("🧬 Analyzing functional content...")
    functional_content = analyze_functional_content(annotations)
    
    print("🚀 Analyzing expansion potential...")
    expansion_potential = analyze_expansion_potential(model, annotations)
    
    print("💡 Generating recommendations...")
    recommendations = generate_expansion_recommendations({
        'model_coverage': model_coverage,
        'functional_content': functional_content,
        'expansion_potential': expansion_potential
    })
    
    # Compile full analysis
    analysis = {
        'metadata': {
            'model_file': args.model,
            'emapper_file': args.emapper,
            'min_score_threshold': args.min_score,
            'max_evalue_threshold': args.max_evalue,
            'total_annotations_analyzed': len(annotations)
        },
        'model_coverage': model_coverage,
        'functional_content': functional_content,
        'expansion_potential': expansion_potential,
        'recommendations': recommendations
    }
    
    # Save analysis
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with output_path.open('w') as f:
        json.dump(analysis, f, indent=2, default=str)
    
    print(f"\n💾 Analysis saved to: {output_path}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("📈 EXPANSION POTENTIAL SUMMARY")
    print("=" * 60)
    
    gene_cov = model_coverage['gene_coverage']
    print(f"Gene Coverage: {gene_cov['coverage_percentage']:.1f}% ({gene_cov['genes_in_common']:,}/{gene_cov['model_genes']:,} genes)")
    print(f"Additional genes in eggNOG: {gene_cov['genes_only_in_emapper']:,}")
    
    func = functional_content
    print(f"\nFunctional Annotations Available:")
    print(f"  • EC numbers: {func['ec_numbers']['unique_ecs']:,} unique ({func['ec_numbers']['genes_with_ec']:,} genes)")
    print(f"  • KEGG reactions: {func['kegg_reactions']['unique_reactions']:,} unique ({func['kegg_reactions']['genes_with_kegg_reaction']:,} genes)")
    print(f"  • KEGG pathways: {func['kegg_pathways']['unique_pathways']:,} unique")
    print(f"  • GO terms: {func['go_terms']['unique_go_terms']:,} unique ({func['go_terms']['genes_with_go']:,} genes)")
    
    gpr_pot = expansion_potential['gpr_completion_potential']
    print(f"\nGPR Completion Potential:")
    print(f"  • Reactions without GPRs: {gpr_pot['reactions_without_gprs']:,}")
    print(f"  • Could complete GPRs for: {gpr_pot['reactions_with_potential_gprs']:,} ({gpr_pot['completion_percentage']:.1f}%)")
    
    print(f"\n🎯 HIGH PRIORITY RECOMMENDATIONS:")
    for i, rec in enumerate(recommendations['high_priority'], 1):
        print(f"  {i}. {rec['action']}")
        print(f"     Impact: {rec['impact']}")
        print(f"     Difficulty: {rec['difficulty']}")
    
    return 0


if __name__ == "__main__":
    exit(main())