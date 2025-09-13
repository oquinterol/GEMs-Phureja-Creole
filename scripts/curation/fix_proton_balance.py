#!/usr/bin/env python3
"""
Fix proton (H+) balance issues in reactions.

Many reactions have H imbalances of ±1 that can be fixed by adding H+ ions.
This script identifies such cases and proposes fixes.

Usage:
  python scripts/fix_proton_balance.py \
    --model models/current/creole_h_balanced.xml \
    --out models/current/creole_proton_balanced.xml
"""
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

import cobra
from cobra.io import read_sbml_model, write_sbml_model


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fix proton balance by adding H+ ions")
    p.add_argument("--model", required=True, help="Input SBML model")
    p.add_argument("--out", required=True, help="Output fixed model")
    p.add_argument("--dry-run", action="store_true", help="Analyze only, don't modify")
    p.add_argument("--verbose", action="store_true", help="Verbose output")
    return p.parse_args()


def parse_formula(formula: str) -> Dict[str, int]:
    """Parse chemical formula into element counts."""
    if not formula or formula == "" or formula.lower() in ['none', 'null']:
        return {}
    
    # Handle charge notation (remove charge for mass balance)
    formula_clean = re.sub(r'[+-]\d*$', '', formula)
    
    # Parse elements using regex
    element_pattern = r'([A-Z][a-z]?)(\d*)'
    elements = {}
    
    for match in re.finditer(element_pattern, formula_clean):
        element = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        elements[element] = elements.get(element, 0) + count
    
    return elements


def analyze_h_balance(reaction: cobra.Reaction) -> Dict:
    """Analyze H balance for a reaction."""
    
    reactant_h = 0
    product_h = 0  
    missing_formulas = []
    
    # Process reactants
    for metabolite in reaction.reactants:
        coefficient = reaction.get_coefficient(metabolite.id)
        if not metabolite.formula:
            missing_formulas.append(metabolite.id)
            continue
            
        elements = parse_formula(metabolite.formula)
        reactant_h += abs(coefficient) * elements.get('H', 0)
    
    # Process products
    for metabolite in reaction.products:
        coefficient = reaction.get_coefficient(metabolite.id)
        if not metabolite.formula:
            missing_formulas.append(metabolite.id)
            continue
            
        elements = parse_formula(metabolite.formula)
        product_h += abs(coefficient) * elements.get('H', 0)
    
    h_imbalance = reactant_h - product_h
    
    return {
        'h_imbalance': h_imbalance,
        'missing_formulas': missing_formulas,
        'can_fix_with_proton': abs(h_imbalance) > 0 and abs(h_imbalance) <= 3,
        'protons_needed': h_imbalance
    }


def find_proton_metabolite(model: cobra.Model) -> Optional[cobra.Metabolite]:
    """Find proton (H+) metabolite in the model."""
    
    # Common proton identifiers
    proton_ids = ['h', 'cpd00067', 'H', 'proton', 'h_plus']
    proton_patterns = [r'.*cpd00067.*', r'^h_.*', r'.*proton.*']
    
    # Direct ID match
    for proton_id in proton_ids:
        for compartment in ['_c0', '_d0', '_m0', '_p0', '_e0', '_x0', '_v0']:
            full_id = proton_id + compartment
            try:
                candidate = model.metabolites.get_by_id(full_id)
                # Verify it's actually a proton by checking formula
                if candidate.formula == 'H' or candidate.formula == 'H+':
                    return candidate
            except KeyError:
                continue
    
    # Pattern match in existing metabolites
    for met in model.metabolites:
        # Check formula first (most reliable)
        if met.formula in ['H', 'H+', 'H1']:
            return met
        
        # Check ID patterns
        for pattern in proton_patterns:
            if re.match(pattern, met.id, re.IGNORECASE):
                return met
        
        # Check name patterns
        if met.name and ('proton' in met.name.lower() or 'h+' in met.name.lower()):
            return met
    
    return None


def fix_proton_balance(model: cobra.Model, dry_run: bool = False, verbose: bool = False) -> int:
    """Fix proton balance issues by adding H+ ions."""
    
    proton_met = find_proton_metabolite(model)
    if not proton_met:
        print("❌ Could not find proton (H+) metabolite in model")
        return 0
    
    print(f"⚡ Found proton metabolite: {proton_met.id} ({proton_met.name})")
    
    fixes_made = 0
    remaining_imbalanced = []
    
    for reaction in model.reactions:
        balance = analyze_h_balance(reaction)
        
        if balance['missing_formulas']:
            continue  # Skip reactions with missing formulas
        
        if balance['can_fix_with_proton'] and abs(balance['protons_needed']) > 0:
            protons_needed = balance['protons_needed']
            
            if verbose:
                print(f"  {reaction.id}: needs {protons_needed} H+ ions")
            
            if not dry_run:
                # Add protons to the appropriate side
                if protons_needed > 0:
                    # Add protons to products (consuming protons on reactant side)
                    reaction.add_metabolites({proton_met: -abs(protons_needed)})
                else:
                    # Add protons to reactants (producing protons on product side) 
                    reaction.add_metabolites({proton_met: abs(protons_needed)})
            
            fixes_made += 1
        
        elif abs(balance['h_imbalance']) > 3:
            # Reactions with large H imbalances that need manual review
            remaining_imbalanced.append({
                'id': reaction.id,
                'h_imbalance': balance['h_imbalance']
            })
    
    print(f"Fixed {fixes_made} reactions with proton balance")
    
    if remaining_imbalanced:
        print(f"Found {len(remaining_imbalanced)} reactions with large H imbalances (>3):")
        for rxn in remaining_imbalanced[:10]:  # Show first 10
            print(f"  {rxn['id']}: H={rxn['h_imbalance']}")
        if len(remaining_imbalanced) > 10:
            print(f"  ... and {len(remaining_imbalanced) - 10} more")
    
    return fixes_made


def main() -> int:
    args = parse_args()
    
    print("⚡ PROTON BALANCE FIX")
    print("=" * 30)
    
    # Load model
    print(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)
    
    print(f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    
    # Fix proton balance
    fixes_made = fix_proton_balance(model, args.dry_run, args.verbose)
    
    if not args.dry_run and fixes_made > 0:
        # Save model
        output_path = Path(args.out)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            write_sbml_model(model, str(output_path))
            print(f"\n💾 Fixed model saved to: {output_path}")
        except Exception as e:
            print(f"❌ Error saving model: {e}")
            return 1
    
    print(f"\n✅ Proton balance fix completed!")
    print(f"Total fixes applied: {fixes_made}")
    
    return 0


if __name__ == "__main__":
    exit(main())