#!/usr/bin/env python3
"""
Analyze and fix stoichiometric inconsistencies in the metabolic model.

This script addresses MEMOTE stoichiometric issues:
1. Analyzes stoichiometric consistency
2. Identifies blocked reactions
3. Finds dead-end metabolites
4. Detects energy-generating cycles
5. Proposes fixes for common issues

Usage:
  python scripts/fix_stoichiometric_issues.py \
    --model models/current/creole_biomass_optimized.xml \
    --out models/current/creole_consistent.xml
"""
from __future__ import annotations

import argparse
import json
import numpy as np
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

import cobra
from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis import find_blocked_reactions
from cobra.util import create_stoichiometric_matrix


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fix stoichiometric inconsistencies")
    p.add_argument("--model", required=True, help="Input SBML model")
    p.add_argument("--out", required=True, help="Output fixed model")
    p.add_argument("--dry-run", action="store_true", help="Analyze only, don't modify")
    p.add_argument("--verbose", action="store_true", help="Verbose output")
    return p.parse_args()


def analyze_stoichiometric_consistency(model: cobra.Model) -> Dict:
    """Analyze stoichiometric matrix consistency."""
    
    print("🔬 Analyzing stoichiometric consistency...")
    
    # Create stoichiometric matrix
    S = create_stoichiometric_matrix(model, array_type='dense')
    
    # Calculate matrix rank
    rank = np.linalg.matrix_rank(S)
    m, n = S.shape  # m = metabolites, n = reactions
    
    # Null space dimension
    null_space_dim = n - rank
    
    print(f"  Matrix dimensions: {m} metabolites × {n} reactions")
    print(f"  Matrix rank: {rank}")
    print(f"  Null space dimension: {null_space_dim}")
    
    # Check for conservation laws (rank deficiency)
    conservation_laws = m - rank
    
    consistency_score = rank / m * 100 if m > 0 else 0
    
    return {
        'matrix_shape': (m, n),
        'rank': rank,
        'null_space_dimension': null_space_dim,
        'conservation_laws': conservation_laws,
        'consistency_score': consistency_score,
        'is_consistent': conservation_laws == 0
    }


def find_deadend_metabolites(model: cobra.Model) -> Dict:
    """Find dead-end metabolites (produced or consumed only)."""
    
    print("🔍 Finding dead-end metabolites...")
    
    producers = defaultdict(set)  # metabolites → reactions that produce them
    consumers = defaultdict(set)  # metabolites → reactions that consume them
    
    for reaction in model.reactions:
        for metabolite, coefficient in reaction.metabolites.items():
            if coefficient > 0:  # Product
                producers[metabolite].add(reaction)
            elif coefficient < 0:  # Reactant
                consumers[metabolite].add(reaction)
    
    # Find dead-ends
    produced_only = []  # Only produced, never consumed
    consumed_only = []  # Only consumed, never produced
    
    for metabolite in model.metabolites:
        is_produced = metabolite in producers
        is_consumed = metabolite in consumers
        
        if is_produced and not is_consumed:
            produced_only.append(metabolite)
        elif is_consumed and not is_produced:
            consumed_only.append(metabolite)
    
    # Exclude exchange/demand/sink metabolites
    exchange_mets = set()
    for reaction in model.reactions:
        if len(reaction.metabolites) == 1:  # Exchange/demand/sink
            exchange_mets.update(reaction.metabolites.keys())
    
    produced_only = [m for m in produced_only if m not in exchange_mets]
    consumed_only = [m for m in consumed_only if m not in exchange_mets]
    
    print(f"  Dead-ends (produced only): {len(produced_only)}")
    print(f"  Dead-ends (consumed only): {len(consumed_only)}")
    
    return {
        'produced_only': [m.id for m in produced_only],
        'consumed_only': [m.id for m in consumed_only],
        'total_deadends': len(produced_only) + len(consumed_only)
    }


def analyze_blocked_reactions(model: cobra.Model) -> Dict:
    """Analyze blocked reactions."""
    
    print("🚫 Analyzing blocked reactions...")
    
    # Find blocked reactions using COBRApy
    blocked = find_blocked_reactions(model)
    
    # Categorize blocked reactions
    categories = {
        'exchange': [],
        'demand': [],
        'sink': [],
        'transport': [],
        'metabolic': []
    }
    
    for rxn_id in blocked:
        rxn = model.reactions.get_by_id(rxn_id)
        
        if rxn_id.startswith('EX_'):
            categories['exchange'].append(rxn_id)
        elif rxn_id.startswith('DM_'):
            categories['demand'].append(rxn_id)
        elif rxn_id.startswith('SK_'):
            categories['sink'].append(rxn_id)
        elif len(set(met.compartment for met in rxn.metabolites)) > 1:
            categories['transport'].append(rxn_id)
        else:
            categories['metabolic'].append(rxn_id)
    
    print(f"  Total blocked reactions: {len(blocked)}")
    for category, rxns in categories.items():
        if rxns:
            print(f"    {category}: {len(rxns)}")
    
    return {
        'blocked_reactions': list(blocked),
        'categories': categories,
        'total_blocked': len(blocked),
        'blocked_percentage': len(blocked) / len(model.reactions) * 100
    }


def detect_energy_cycles(model: cobra.Model) -> Dict:
    """Detect potential energy-generating cycles."""
    
    print("⚡ Detecting energy-generating cycles...")
    
    # Focus on energy-carrying metabolites
    energy_metabolites = [
        'cpd00002',  # ATP
        'cpd00008',  # ADP
        'cpd00006',  # NADP+
        'cpd00005',  # NADPH
        'cpd00004',  # NADH
        'cpd00003',  # NAD+
    ]
    
    energy_imbalances = {}
    
    for met_pattern in energy_metabolites:
        # Find all compartments for this metabolite
        mets = [met for met in model.metabolites if met_pattern in met.id]
        
        for met in mets:
            # Count reactions producing vs consuming
            producing = sum(1 for rxn in met.reactions if rxn.metabolites[met] > 0)
            consuming = sum(1 for rxn in met.reactions if rxn.metabolites[met] < 0)
            
            if producing > 0 and consuming > 0:
                ratio = producing / consuming
                energy_imbalances[met.id] = {
                    'producing_reactions': producing,
                    'consuming_reactions': consuming,
                    'ratio': ratio
                }
    
    # Flag potentially problematic ratios
    suspicious = {met_id: info for met_id, info in energy_imbalances.items() 
                 if info['ratio'] > 2.0 or info['ratio'] < 0.5}
    
    print(f"  Energy metabolites analyzed: {len(energy_imbalances)}")
    print(f"  Suspicious ratios: {len(suspicious)}")
    
    return {
        'energy_imbalances': energy_imbalances,
        'suspicious_metabolites': suspicious
    }


def fix_common_issues(model: cobra.Model, analysis_results: Dict, dry_run: bool = False) -> int:
    """Fix common stoichiometric issues."""
    
    print("🔧 Attempting to fix common issues...")
    
    fixes_applied = 0
    
    # 1. Remove duplicate reactions (if any)
    print("  Checking for duplicate reactions...")
    reaction_strings = {}
    duplicates = []
    
    for rxn in model.reactions:
        rxn_string = rxn.build_reaction_string()
        if rxn_string in reaction_strings:
            duplicates.append(rxn)
        else:
            reaction_strings[rxn_string] = rxn
    
    if duplicates and not dry_run:
        print(f"    Removing {len(duplicates)} duplicate reactions")
        model.remove_reactions(duplicates)
        fixes_applied += len(duplicates)
    else:
        print(f"    Found {len(duplicates)} duplicate reactions")
    
    # 2. Add missing reversibility constraints
    print("  Checking reaction reversibility...")
    reversibility_fixes = 0
    
    for rxn in model.reactions:
        # If lower_bound is 0 but reaction involves only reversible biochemistry
        if rxn.lower_bound == 0 and rxn.upper_bound > 0:
            # Check if this should be reversible based on reaction type
            if not rxn.id.startswith(('EX_', 'DM_', 'SK_')):
                # Most metabolic reactions should be reversible unless specified
                if not dry_run:
                    rxn.lower_bound = -1000.0
                    reversibility_fixes += 1
    
    if reversibility_fixes > 0:
        print(f"    Fixed reversibility for {reversibility_fixes} reactions")
        fixes_applied += reversibility_fixes
    
    # 3. Fix obvious transport reactions
    print("  Checking transport reactions...")
    transport_fixes = 0
    
    for rxn in model.reactions:
        if len(set(met.compartment for met in rxn.metabolites)) > 1:
            # This is a transport reaction
            if rxn.lower_bound == 0 and rxn.upper_bound == 1000:
                # Make bidirectional transport
                if not dry_run:
                    rxn.lower_bound = -1000.0
                    transport_fixes += 1
    
    if transport_fixes > 0:
        print(f"    Fixed transport directionality for {transport_fixes} reactions")
        fixes_applied += transport_fixes
    
    return fixes_applied


def generate_consistency_report(analysis_results: Dict, output_path: Path):
    """Generate stoichiometric consistency report."""
    
    report_path = output_path.parent / f"{output_path.stem}_consistency_report.json"
    
    with report_path.open('w') as f:
        json.dump(analysis_results, f, indent=2, default=str)
    
    print(f"📊 Consistency report saved to: {report_path}")
    
    # Generate human-readable summary
    summary_path = output_path.parent / f"{output_path.stem}_consistency_summary.txt"
    
    with summary_path.open('w') as f:
        f.write("STOICHIOMETRIC CONSISTENCY ANALYSIS\\n")
        f.write("=" * 50 + "\\n\\n")
        
        # Stoichiometric matrix
        if 'stoichiometric' in analysis_results:
            stoi = analysis_results['stoichiometric']
            f.write(f"Stoichiometric Matrix:\\n")
            f.write(f"  Shape: {stoi['matrix_shape'][0]} metabolites × {stoi['matrix_shape'][1]} reactions\\n")
            f.write(f"  Rank: {stoi['rank']}\\n")
            f.write(f"  Consistency score: {stoi['consistency_score']:.1f}%\\n")
            f.write(f"  Conservation laws: {stoi['conservation_laws']}\\n\\n")
        
        # Blocked reactions
        if 'blocked' in analysis_results:
            blocked = analysis_results['blocked']
            f.write(f"Blocked Reactions:\\n")
            f.write(f"  Total blocked: {blocked['total_blocked']} ({blocked['blocked_percentage']:.1f}%)\\n")
            for category, count in blocked['categories'].items():
                if count:
                    f.write(f"    {category}: {len(count)}\\n")
            f.write("\\n")
        
        # Dead-ends
        if 'deadends' in analysis_results:
            dead = analysis_results['deadends']
            f.write(f"Dead-end Metabolites:\\n")
            f.write(f"  Produced only: {len(dead['produced_only'])}\\n")
            f.write(f"  Consumed only: {len(dead['consumed_only'])}\\n")
            f.write(f"  Total dead-ends: {dead['total_deadends']}\\n\\n")
        
        # Energy cycles
        if 'energy_cycles' in analysis_results:
            energy = analysis_results['energy_cycles']
            f.write(f"Energy Cycle Analysis:\\n")
            f.write(f"  Suspicious metabolites: {len(energy['suspicious_metabolites'])}\\n")
            for met_id, info in list(energy['suspicious_metabolites'].items())[:5]:
                f.write(f"    {met_id}: ratio {info['ratio']:.2f}\\n")
    
    print(f"📝 Consistency summary saved to: {summary_path}")


def main() -> int:
    args = parse_args()
    
    print("⚖️ STOICHIOMETRIC CONSISTENCY ANALYSIS AND FIXES")
    print("=" * 60)
    
    # Load model
    print(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)
    
    print(f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    
    # Perform analyses
    analysis_results = {}
    
    # 1. Stoichiometric consistency
    analysis_results['stoichiometric'] = analyze_stoichiometric_consistency(model)
    
    # 2. Dead-end metabolites
    analysis_results['deadends'] = find_deadend_metabolites(model)
    
    # 3. Blocked reactions
    analysis_results['blocked'] = analyze_blocked_reactions(model)
    
    # 4. Energy cycles
    analysis_results['energy_cycles'] = detect_energy_cycles(model)
    
    # 5. Apply fixes
    fixes_applied = fix_common_issues(model, analysis_results, args.dry_run)
    analysis_results['fixes_applied'] = fixes_applied
    
    # Print summary
    print("\\n" + "=" * 60)
    print("📊 CONSISTENCY ANALYSIS SUMMARY")
    print("=" * 60)
    
    stoi = analysis_results['stoichiometric']
    print(f"Stoichiometric consistency: {stoi['consistency_score']:.1f}%")
    print(f"Conservation laws: {stoi['conservation_laws']}")
    
    deadends = analysis_results['deadends']
    print(f"Dead-end metabolites: {deadends['total_deadends']}")
    
    blocked = analysis_results['blocked']
    print(f"Blocked reactions: {blocked['total_blocked']} ({blocked['blocked_percentage']:.1f}%)")
    
    energy = analysis_results['energy_cycles']
    print(f"Suspicious energy ratios: {len(energy['suspicious_metabolites'])}")
    
    print(f"Fixes applied: {fixes_applied}")
    
    # Save results
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if not args.dry_run and fixes_applied > 0:
        try:
            write_sbml_model(model, str(output_path))
            print(f"\\n💾 Fixed model saved to: {output_path}")
        except Exception as e:
            print(f"❌ Error saving model: {e}")
            return 1
    
    # Generate report
    generate_consistency_report(analysis_results, output_path)
    
    print(f"\\n✅ Stoichiometric analysis completed!")
    
    return 0


if __name__ == "__main__":
    exit(main())