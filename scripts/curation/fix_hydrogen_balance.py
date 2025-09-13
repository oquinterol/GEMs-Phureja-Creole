#!/usr/bin/env python3
"""
Fix hydrogen balance issues by adding water molecules.

Many reactions are missing H2O molecules to balance H and O atoms.
This script identifies such cases and proposes fixes.

Usage:
  python scripts/fix_hydrogen_balance.py \
    --model models/current/creole_balanced.xml \
    --out models/current/creole_h_balanced.xml
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
    p = argparse.ArgumentParser(description="Fix hydrogen balance by adding water molecules")
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


def analyze_h_o_balance(reaction: cobra.Reaction) -> Dict:
    """Analyze H and O balance for a reaction."""
    
    reactant_h = 0
    reactant_o = 0
    product_h = 0  
    product_o = 0
    missing_formulas = []
    
    # Process reactants
    for metabolite in reaction.reactants:
        coefficient = reaction.get_coefficient(metabolite.id)
        if not metabolite.formula:
            missing_formulas.append(metabolite.id)
            continue
            
        elements = parse_formula(metabolite.formula)
        reactant_h += abs(coefficient) * elements.get('H', 0)
        reactant_o += abs(coefficient) * elements.get('O', 0)
    
    # Process products
    for metabolite in reaction.products:
        coefficient = reaction.get_coefficient(metabolite.id)
        if not metabolite.formula:
            missing_formulas.append(metabolite.id)
            continue
            
        elements = parse_formula(metabolite.formula)
        product_h += abs(coefficient) * elements.get('H', 0)
        product_o += abs(coefficient) * elements.get('O', 0)
    
    h_imbalance = reactant_h - product_h
    o_imbalance = reactant_o - product_o
    
    return {
        'h_imbalance': h_imbalance,
        'o_imbalance': o_imbalance,
        'missing_formulas': missing_formulas,
        'can_fix_with_water': h_imbalance == 2 * o_imbalance and abs(h_imbalance) > 0,
        'water_molecules_needed': o_imbalance if h_imbalance == 2 * o_imbalance else 0
    }


def find_water_metabolite(model: cobra.Model) -> Optional[cobra.Metabolite]:
    """Find water metabolite in the model."""
    
    # Common water identifiers
    water_ids = ['h2o', 'cpd00001', 'H2O', 'water']
    water_patterns = [r'.*h2o.*', r'.*water.*', r'.*cpd00001.*']
    
    # Direct ID match
    for water_id in water_ids:
        for compartment in ['_c0', '_d0', '_m0', '_p0', '_e0', '_x0', '_v0']:
            full_id = water_id + compartment
            try:
                return model.metabolites.get_by_id(full_id)
            except KeyError:
                continue
    
    # Pattern match in existing metabolites
    for met in model.metabolites:
        # Check formula first (most reliable)
        if met.formula == 'H2O':
            return met
        
        # Check ID patterns
        for pattern in water_patterns:
            if re.match(pattern, met.id, re.IGNORECASE):
                return met
        
        # Check name patterns
        if met.name and 'water' in met.name.lower():
            return met
    
    return None


def fix_hydrogen_balance_with_water(model: cobra.Model, dry_run: bool = False, verbose: bool = False) -> int:
    """Fix hydrogen balance issues by adding water molecules."""
    
    water_met = find_water_metabolite(model)
    if not water_met:
        print("❌ Could not find water metabolite in model")
        return 0
    
    print(f"🌊 Found water metabolite: {water_met.id} ({water_met.name})")
    
    fixes_made = 0
    problematic_reactions = []
    
    for reaction in model.reactions:
        balance = analyze_h_o_balance(reaction)
        
        if balance['missing_formulas']:
            continue  # Skip reactions with missing formulas
        
        if balance['can_fix_with_water'] and abs(balance['water_molecules_needed']) > 0:
            water_needed = balance['water_molecules_needed']
            
            if verbose:
                print(f"  {reaction.id}: needs {water_needed} H2O molecules")
            
            if not dry_run:
                # Add water to the appropriate side
                if water_needed > 0:
                    # Add water to products (consuming water on reactant side)
                    reaction.add_metabolites({water_met: -abs(water_needed)})
                else:
                    # Add water to reactants (producing water on product side) 
                    reaction.add_metabolites({water_met: abs(water_needed)})
            
            fixes_made += 1
        
        elif abs(balance['h_imbalance']) > 0 or abs(balance['o_imbalance']) > 0:
            # Reactions that can't be fixed with just water
            problematic_reactions.append({
                'id': reaction.id,
                'h_imbalance': balance['h_imbalance'],
                'o_imbalance': balance['o_imbalance']
            })
    
    print(f"Fixed {fixes_made} reactions with water molecules")
    
    if problematic_reactions:
        print(f"Found {len(problematic_reactions)} reactions that need manual review:")
        for rxn in problematic_reactions[:10]:  # Show first 10
            print(f"  {rxn['id']}: H={rxn['h_imbalance']}, O={rxn['o_imbalance']}")
        if len(problematic_reactions) > 10:
            print(f"  ... and {len(problematic_reactions) - 10} more")
    
    return fixes_made


def main() -> int:
    args = parse_args()
    
    print("🌊 HYDROGEN BALANCE FIX WITH WATER")
    print("=" * 50)
    
    # Load model
    print(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)
    
    print(f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    
    # Fix hydrogen balance
    fixes_made = fix_hydrogen_balance_with_water(model, args.dry_run, args.verbose)
    
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
    
    print(f"\n✅ Hydrogen balance fix completed!")
    print(f"Total fixes applied: {fixes_made}")
    
    return 0


if __name__ == "__main__":
    exit(main())