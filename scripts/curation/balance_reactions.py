#!/usr/bin/env python3
"""
Balance reactions for mass and charge consistency.

This script analyzes and fixes mass/charge balance issues in the model:
1. Identifies metabolites without formulas
2. Attempts to fill missing formulas using ChEBI mappings
3. Analyzes mass and charge balance for all reactions
4. Reports imbalanced reactions
5. Provides suggestions for fixes

Usage:
  python scripts/balance_reactions.py \
    --model models/current/creole_unified.xml \
    --out models/current/creole_balanced.xml
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
import requests
import time


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Balance reactions for mass/charge consistency")
    p.add_argument("--model", required=True, help="Input SBML model")
    p.add_argument("--out", required=True, help="Output balanced model")
    p.add_argument("--verbose", action="store_true", help="Verbose output")
    p.add_argument("--dry-run", action="store_true", help="Analyze only, don't modify")
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


def get_charge_from_formula(formula: str) -> int:
    """Extract charge from formula."""
    if not formula:
        return 0
    
    charge_match = re.search(r'([+-])(\d*)$', formula)
    if charge_match:
        sign = 1 if charge_match.group(1) == '+' else -1
        magnitude = int(charge_match.group(2)) if charge_match.group(2) else 1
        return sign * magnitude
    
    return 0


def fetch_chebi_formula(chebi_id: str) -> Optional[Tuple[str, int]]:
    """Fetch formula and charge from ChEBI API."""
    try:
        # Rate limiting
        time.sleep(0.1)
        
        url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{chebi_id}"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            # This is a simplified parser - in production would use proper XML/HTML parsing
            text = response.text
            
            # Look for molecular formula
            formula_match = re.search(r'Molecular Formula</td[^>]*>\s*<td[^>]*>([^<]+)</td>', text)
            formula = formula_match.group(1).strip() if formula_match else None
            
            # Look for charge
            charge_match = re.search(r'Charge</td[^>]*>\s*<td[^>]*>([^<]+)</td>', text)
            charge = 0
            if charge_match:
                try:
                    charge = int(charge_match.group(1).strip())
                except ValueError:
                    pass
            
            return formula, charge
            
    except Exception as e:
        print(f"  Warning: Failed to fetch ChEBI:{chebi_id}: {e}")
        return None


def enrich_metabolite_formulas(model: cobra.Model, verbose: bool = False) -> int:
    """Enrich missing metabolite formulas using ChEBI annotations."""
    
    formulas_added = 0
    
    for metabolite in model.metabolites:
        if metabolite.formula:
            continue  # Skip if formula already present
        
        # Look for ChEBI annotation
        chebi_ids = []
        if hasattr(metabolite, 'annotation'):
            for key, value in metabolite.annotation.items():
                if 'chebi' in key.lower():
                    if isinstance(value, list):
                        for v in value:
                            if isinstance(v, str):
                                # Extract ChEBI ID
                                chebi_match = re.search(r'(?:CHEBI:)?(\d+)', v)
                                if chebi_match:
                                    chebi_ids.append(chebi_match.group(1))
                    elif isinstance(value, str):
                        chebi_match = re.search(r'(?:CHEBI:)?(\d+)', value)
                        if chebi_match:
                            chebi_ids.append(chebi_match.group(1))
        
        # Try to fetch formula from ChEBI
        for chebi_id in chebi_ids[:1]:  # Try first ChEBI ID only
            if verbose:
                print(f"  Fetching formula for {metabolite.id} from ChEBI:{chebi_id}")
            
            result = fetch_chebi_formula(chebi_id)
            if result:
                formula, charge = result
                if formula:
                    metabolite.formula = formula
                    if charge != 0 and not hasattr(metabolite, 'charge'):
                        metabolite.charge = charge
                    formulas_added += 1
                    if verbose:
                        print(f"    Added formula: {formula} (charge: {charge})")
                    break
    
    return formulas_added


def analyze_reaction_balance(reaction: cobra.Reaction) -> Dict:
    """Analyze mass and charge balance for a reaction."""
    
    # Calculate element balance
    reactant_elements = defaultdict(float)
    product_elements = defaultdict(float)
    reactant_charge = 0
    product_charge = 0
    
    missing_formulas = []
    
    # Process reactants
    for metabolite in reaction.reactants:
        coefficient = reaction.get_coefficient(metabolite.id)
        if not metabolite.formula:
            missing_formulas.append(metabolite.id)
            continue
            
        elements = parse_formula(metabolite.formula)
        charge = get_charge_from_formula(metabolite.formula)
        
        for element, count in elements.items():
            reactant_elements[element] += abs(coefficient) * count
        
        reactant_charge += abs(coefficient) * charge
    
    # Process products
    for metabolite in reaction.products:
        coefficient = reaction.get_coefficient(metabolite.id)
        if not metabolite.formula:
            missing_formulas.append(metabolite.id)
            continue
            
        elements = parse_formula(metabolite.formula)
        charge = get_charge_from_formula(metabolite.formula)
        
        for element, count in elements.items():
            product_elements[element] += abs(coefficient) * count
        
        product_charge += abs(coefficient) * charge
    
    # Calculate imbalances
    all_elements = set(reactant_elements.keys()) | set(product_elements.keys())
    element_imbalances = {}
    
    for element in all_elements:
        reactant_count = reactant_elements.get(element, 0)
        product_count = product_elements.get(element, 0)
        imbalance = reactant_count - product_count
        if abs(imbalance) > 1e-6:
            element_imbalances[element] = imbalance
    
    charge_imbalance = reactant_charge - product_charge
    
    return {
        'missing_formulas': missing_formulas,
        'element_imbalances': element_imbalances,
        'charge_imbalance': charge_imbalance,
        'is_mass_balanced': len(element_imbalances) == 0,
        'is_charge_balanced': abs(charge_imbalance) < 1e-6,
        'is_balanced': len(element_imbalances) == 0 and abs(charge_imbalance) < 1e-6
    }


def analyze_model_balance(model: cobra.Model, verbose: bool = False) -> Dict:
    """Analyze balance for all reactions in the model."""
    
    print("🔬 Analyzing reaction balance...")
    
    results = {
        'total_reactions': len(model.reactions),
        'reactions_with_missing_formulas': 0,
        'mass_imbalanced_reactions': 0,
        'charge_imbalanced_reactions': 0,
        'fully_balanced_reactions': 0,
        'reaction_details': {},
        'missing_formula_metabolites': set(),
        'element_imbalance_summary': Counter(),
        'charge_imbalance_distribution': []
    }
    
    for reaction in model.reactions:
        balance = analyze_reaction_balance(reaction)
        results['reaction_details'][reaction.id] = balance
        
        # Update statistics
        if balance['missing_formulas']:
            results['reactions_with_missing_formulas'] += 1
            results['missing_formula_metabolites'].update(balance['missing_formulas'])
        
        if not balance['is_mass_balanced']:
            results['mass_imbalanced_reactions'] += 1
            for element, imbalance in balance['element_imbalances'].items():
                results['element_imbalance_summary'][element] += 1
        
        if not balance['is_charge_balanced']:
            results['charge_imbalanced_reactions'] += 1
            results['charge_imbalance_distribution'].append(balance['charge_imbalance'])
        
        if balance['is_balanced']:
            results['fully_balanced_reactions'] += 1
        
        if verbose and not balance['is_balanced']:
            print(f"  Imbalanced: {reaction.id}")
            if balance['element_imbalances']:
                print(f"    Elements: {balance['element_imbalances']}")
            if not balance['is_charge_balanced']:
                print(f"    Charge: {balance['charge_imbalance']}")
    
    return results


def generate_balance_report(results: Dict, output_path: Path):
    """Generate detailed balance analysis report."""
    
    report_path = output_path.parent / f"{output_path.stem}_balance_report.json"
    
    # Convert sets to lists for JSON serialization
    report_data = results.copy()
    report_data['missing_formula_metabolites'] = list(results['missing_formula_metabolites'])
    
    with report_path.open('w') as f:
        json.dump(report_data, f, indent=2, default=str)
    
    print(f"📊 Balance report saved to: {report_path}")
    
    # Generate human-readable summary
    summary_path = output_path.parent / f"{output_path.stem}_balance_summary.txt"
    
    with summary_path.open('w') as f:
        f.write("REACTION BALANCE ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Total reactions analyzed: {results['total_reactions']:,}\n")
        f.write(f"Fully balanced reactions: {results['fully_balanced_reactions']:,} ({results['fully_balanced_reactions']/results['total_reactions']*100:.1f}%)\n")
        f.write(f"Mass imbalanced reactions: {results['mass_imbalanced_reactions']:,} ({results['mass_imbalanced_reactions']/results['total_reactions']*100:.1f}%)\n")
        f.write(f"Charge imbalanced reactions: {results['charge_imbalanced_reactions']:,} ({results['charge_imbalanced_reactions']/results['total_reactions']*100:.1f}%)\n")
        f.write(f"Reactions with missing formulas: {results['reactions_with_missing_formulas']:,} ({results['reactions_with_missing_formulas']/results['total_reactions']*100:.1f}%)\n\n")
        
        f.write(f"Metabolites without formulas: {len(results['missing_formula_metabolites']):,}\n\n")
        
        f.write("Most common element imbalances:\n")
        for element, count in results['element_imbalance_summary'].most_common(10):
            f.write(f"  {element}: {count} reactions\n")
        
        f.write("\nTop 20 most problematic reactions:\n")
        problematic = []
        for rxn_id, details in results['reaction_details'].items():
            if not details['is_balanced']:
                score = len(details['element_imbalances']) + (1 if not details['is_charge_balanced'] else 0)
                problematic.append((rxn_id, score, details))
        
        problematic.sort(key=lambda x: x[1], reverse=True)
        for rxn_id, score, details in problematic[:20]:
            f.write(f"  {rxn_id}: {len(details['element_imbalances'])} element imbalances")
            if not details['is_charge_balanced']:
                f.write(f", charge imbalance: {details['charge_imbalance']}")
            f.write("\n")
    
    print(f"📝 Balance summary saved to: {summary_path}")


def main() -> int:
    args = parse_args()
    
    print("⚖️ REACTION BALANCE ANALYSIS AND CORRECTION")
    print("=" * 60)
    
    # Load model
    print(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)
    
    print(f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    
    # Count metabolites without formulas
    missing_formulas = [m for m in model.metabolites if not m.formula]
    print(f"Metabolites without formulas: {len(missing_formulas)}")
    
    if not args.dry_run:
        # Try to enrich formulas using ChEBI
        print("\n🧪 Enriching metabolite formulas from ChEBI...")
        formulas_added = enrich_metabolite_formulas(model, args.verbose)
        print(f"Added {formulas_added} formulas from ChEBI")
        
        # Recount after enrichment
        missing_formulas_after = [m for m in model.metabolites if not m.formula]
        print(f"Metabolites still without formulas: {len(missing_formulas_after)}")
    
    # Analyze balance
    print("\n⚖️ Analyzing reaction balance...")
    balance_results = analyze_model_balance(model, args.verbose)
    
    # Print summary
    print("\n" + "=" * 60)
    print("📊 BALANCE ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"Total reactions: {balance_results['total_reactions']:,}")
    print(f"Fully balanced: {balance_results['fully_balanced_reactions']:,} ({balance_results['fully_balanced_reactions']/balance_results['total_reactions']*100:.1f}%)")
    print(f"Mass imbalanced: {balance_results['mass_imbalanced_reactions']:,} ({balance_results['mass_imbalanced_reactions']/balance_results['total_reactions']*100:.1f}%)")
    print(f"Charge imbalanced: {balance_results['charge_imbalanced_reactions']:,} ({balance_results['charge_imbalanced_reactions']/balance_results['total_reactions']*100:.1f}%)")
    print(f"Missing formulas: {balance_results['reactions_with_missing_formulas']:,} reactions affected")
    print(f"Unique metabolites without formulas: {len(balance_results['missing_formula_metabolites'])}")
    
    print(f"\nMost problematic elements:")
    for element, count in balance_results['element_imbalance_summary'].most_common(5):
        print(f"  {element}: {count} reactions imbalanced")
    
    # Save results
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Generate reports
    generate_balance_report(balance_results, output_path)
    
    if not args.dry_run:
        # Save model (even if no major fixes, formulas might have been added)
        try:
            write_sbml_model(model, str(output_path))
            print(f"\n💾 Model saved to: {output_path}")
        except Exception as e:
            print(f"❌ Error saving model: {e}")
            return 1
    
    print(f"\n✅ Balance analysis completed!")
    
    # Provide recommendations
    print(f"\n💡 RECOMMENDATIONS:")
    if balance_results['missing_formula_metabolites']:
        print(f"1. Add formulas for {len(balance_results['missing_formula_metabolites'])} metabolites")
    if balance_results['mass_imbalanced_reactions'] > 0:
        print(f"2. Review {balance_results['mass_imbalanced_reactions']} mass-imbalanced reactions")
    if balance_results['charge_imbalanced_reactions'] > 0:
        print(f"3. Review {balance_results['charge_imbalanced_reactions']} charge-imbalanced reactions")
    
    return 0


if __name__ == "__main__":
    exit(main())