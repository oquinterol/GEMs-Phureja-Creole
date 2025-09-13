#!/usr/bin/env python3
"""
Optimize biomass reactions and add NGAM (Non-Growth Associated Maintenance).

This script addresses MEMOTE biomass and NGAM issues:
1. Reduces excessive GAM (Growth Associated Maintenance) in biomass
2. Adds proper NGAM reaction
3. Balances biomass composition
4. Validates biomass functionality

Usage:
  python scripts/optimize_biomass_ngam.py \
    --model models/current/creole_proton_balanced.xml \
    --out models/current/creole_biomass_optimized.xml
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
from cobra import Reaction, Metabolite


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Optimize biomass and add NGAM")
    p.add_argument("--model", required=True, help="Input SBML model")
    p.add_argument("--out", required=True, help="Output optimized model")
    p.add_argument("--gam-atp", type=float, default=10.0, help="ATP cost for GAM (default: 10)")
    p.add_argument("--ngam-flux", type=float, default=1.0, help="NGAM maintenance flux (default: 1.0)")
    p.add_argument("--dry-run", action="store_true", help="Analyze only, don't modify")
    p.add_argument("--verbose", action="store_true", help="Verbose output")
    return p.parse_args()


def find_metabolites_by_id_pattern(model: cobra.Model, patterns: List[str]) -> Dict[str, cobra.Metabolite]:
    """Find metabolites by ID pattern."""
    found = {}
    
    for pattern in patterns:
        for met in model.metabolites:
            if pattern in met.id.lower():
                found[pattern] = met
                break
    
    return found


def create_ngam_reaction(model: cobra.Model, flux: float = 1.0) -> Optional[cobra.Reaction]:
    """Create NGAM (Non-Growth Associated Maintenance) reaction."""
    
    # Find required metabolites for ATP hydrolysis
    metabolites_needed = {
        'atp': ['cpd00002_c0', 'atp_c0'],
        'adp': ['cpd00008_c0', 'adp_c0'], 
        'pi': ['cpd00009_c0', 'pi_c0', 'phosphate_c0'],
        'h2o': ['cpd00001_c0', 'h2o_c0'],
        'h': ['cpd00067_c0', 'h_c0']
    }
    
    metabolites = {}
    
    for met_type, candidates in metabolites_needed.items():
        found = None
        for candidate in candidates:
            try:
                found = model.metabolites.get_by_id(candidate)
                break
            except KeyError:
                continue
        
        if found:
            metabolites[met_type] = found
        else:
            print(f"❌ Could not find {met_type} metabolite")
            return None
    
    # Create NGAM reaction: ATP + H2O → ADP + Pi + H+
    ngam = Reaction('NGAM')
    ngam.name = 'Non-growth associated maintenance'
    ngam.subsystem = 'Maintenance'
    ngam.lower_bound = flux  # Fixed flux for maintenance
    ngam.upper_bound = flux
    
    # Add metabolites with stoichiometry
    ngam.add_metabolites({
        metabolites['atp']: -1.0,    # ATP consumed
        metabolites['h2o']: -1.0,    # H2O consumed
        metabolites['adp']: 1.0,     # ADP produced
        metabolites['pi']: 1.0,      # Pi produced
        metabolites['h']: 1.0        # H+ produced
    })
    
    # Add SBO annotation
    ngam.annotation = {'sbo': 'SBO:0000630'}  # ATP maintenance
    
    print(f"🔧 Created NGAM reaction: {ngam.reaction}")
    return ngam


def optimize_biomass_gam(model: cobra.Model, biomass_id: str, target_gam: float = 10.0) -> bool:
    """Optimize GAM (Growth Associated Maintenance) in biomass reaction."""
    
    try:
        biomass = model.reactions.get_by_id(biomass_id)
    except KeyError:
        print(f"❌ Biomass reaction {biomass_id} not found")
        return False
    
    print(f"🎯 Optimizing GAM in {biomass_id}")
    
    # Find ATP/ADP/Pi in biomass
    atp_coef = 0
    adp_coef = 0
    pi_coef = 0
    h_coef = 0
    
    metabolite_coeffs = dict(biomass.metabolites)
    
    for met, coef in metabolite_coeffs.items():
        met_id_lower = met.id.lower()
        if 'cpd00002' in met_id_lower or 'atp' in met_id_lower:  # ATP
            atp_coef = coef
        elif 'cpd00008' in met_id_lower or 'adp' in met_id_lower:  # ADP
            adp_coef = coef
        elif 'cpd00009' in met_id_lower or 'pi' in met_id_lower:  # Pi
            pi_coef = coef
        elif 'cpd00067' in met_id_lower and '_c0' in met.id:  # H+ in cytoplasm
            h_coef = coef
    
    current_gam = abs(atp_coef)
    print(f"  Current GAM: {current_gam} ATP")
    print(f"  Target GAM: {target_gam} ATP")
    
    if current_gam == 0:
        print("  ⚠️ No ATP found in biomass - GAM not present")
        return False
    
    # Calculate adjustment factor
    adjustment_factor = target_gam / current_gam
    
    print(f"  Adjustment factor: {adjustment_factor:.3f}")
    
    # Update ATP/ADP/Pi/H+ coefficients
    updates = {}
    for met, coef in metabolite_coeffs.items():
        met_id_lower = met.id.lower()
        if ('cpd00002' in met_id_lower or 'atp' in met_id_lower or 
            'cpd00008' in met_id_lower or 'adp' in met_id_lower or
            'cpd00009' in met_id_lower or 'pi' in met_id_lower or
            ('cpd00067' in met_id_lower and '_c0' in met.id)):
            
            new_coef = coef * adjustment_factor
            updates[met] = new_coef - coef  # Difference to add
    
    # Apply updates
    biomass.add_metabolites(updates)
    
    print(f"  ✅ Updated GAM from {current_gam} to {target_gam} ATP")
    return True


def validate_biomass_functionality(model: cobra.Model, biomass_id: str) -> Dict:
    """Validate that biomass reaction can produce growth."""
    
    try:
        biomass = model.reactions.get_by_id(biomass_id)
    except KeyError:
        return {'functional': False, 'error': f'Biomass {biomass_id} not found'}
    
    # Test FBA with biomass objective
    original_objective = model.objective
    model.objective = biomass_id
    
    try:
        solution = model.optimize()
        
        result = {
            'functional': solution.status == 'optimal',
            'objective_value': solution.objective_value if solution.status == 'optimal' else 0,
            'status': solution.status
        }
        
        # Check for major precursors
        precursor_check = check_essential_precursors(biomass)
        result.update(precursor_check)
        
    except Exception as e:
        result = {'functional': False, 'error': str(e)}
    finally:
        model.objective = original_objective
    
    return result


def check_essential_precursors(biomass_reaction: cobra.Reaction) -> Dict:
    """Check for essential biomass precursors."""
    
    essential_compounds = {
        'amino_acids': ['cpd00023', 'cpd00033', 'cpd00035', 'cpd00039', 'cpd00041', 'cpd00051', 
                       'cpd00053', 'cpd00054', 'cpd00060', 'cpd00062', 'cpd00065', 'cpd00066',
                       'cpd00069', 'cpd00084', 'cpd00107', 'cpd00119', 'cpd00129', 'cpd00132',
                       'cpd00156', 'cpd00159', 'cpd00161'],
        'nucleotides': ['cpd00002', 'cpd00008', 'cpd00038', 'cpd00052'],
        'lipids': ['cpd00205', 'cpd00214'],
        'carbohydrates': ['cpd00163', 'cpd00099'],
        'cofactors': ['cpd00003', 'cpd00004', 'cpd00005', 'cpd00006', 'cpd00010']
    }
    
    found_categories = {}
    biomass_metabolites = {met.id.split('_')[0] for met in biomass_reaction.reactants}
    
    for category, compounds in essential_compounds.items():
        found = sum(1 for compound in compounds if compound in biomass_metabolites)
        total = len(compounds)
        found_categories[category] = {'found': found, 'total': total, 'percent': found/total*100}
    
    return {'precursor_coverage': found_categories}


def generate_biomass_report(model: cobra.Model, output_path: Path) -> Dict:
    """Generate biomass optimization report."""
    
    biomass_reactions = [rxn for rxn in model.reactions if 'biomass' in rxn.id.lower() or 'bio' in rxn.id.lower()]
    
    report = {
        'biomass_reactions': [],
        'ngam_present': any('ngam' in rxn.id.lower() for rxn in model.reactions),
        'objective_reaction': str(model.objective.expression)
    }
    
    for biomass in biomass_reactions:
        validation = validate_biomass_functionality(model, biomass.id)
        
        biomass_info = {
            'id': biomass.id,
            'name': biomass.name,
            'n_metabolites': len(biomass.metabolites),
            'bounds': [biomass.lower_bound, biomass.upper_bound],
            'validation': validation
        }
        
        report['biomass_reactions'].append(biomass_info)
    
    # Save report
    report_path = output_path.parent / f"{output_path.stem}_biomass_report.json"
    with report_path.open('w') as f:
        json.dump(report, f, indent=2, default=str)
    
    print(f"📊 Biomass report saved to: {report_path}")
    return report


def main() -> int:
    args = parse_args()
    
    print("🧬 BIOMASS AND NGAM OPTIMIZATION")
    print("=" * 50)
    
    # Load model
    print(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)
    
    print(f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    
    # Check current status
    biomass_reactions = [rxn for rxn in model.reactions if 'biomass' in rxn.id.lower() or 'bio' in rxn.id.lower()]
    ngam_reactions = [rxn for rxn in model.reactions if 'ngam' in rxn.id.lower()]
    
    print(f"Found {len(biomass_reactions)} biomass reactions")
    print(f"Found {len(ngam_reactions)} NGAM reactions")
    
    changes_made = 0
    
    if not args.dry_run:
        # Add NGAM if missing
        if not ngam_reactions:
            print("\\n🔧 Adding NGAM reaction...")
            ngam = create_ngam_reaction(model, args.ngam_flux)
            if ngam:
                model.add_reactions([ngam])
                changes_made += 1
                print(f"✅ Added NGAM with flux {args.ngam_flux}")
            else:
                print("❌ Failed to create NGAM reaction")
        
        # Optimize GAM in biomass reactions
        for biomass in biomass_reactions:
            print(f"\\n🎯 Optimizing GAM in {biomass.id}...")
            if optimize_biomass_gam(model, biomass.id, args.gam_atp):
                changes_made += 1
    
    # Validate biomass functionality
    print("\\n🧪 Validating biomass functionality...")
    for biomass in biomass_reactions:
        validation = validate_biomass_functionality(model, biomass.id)
        
        print(f"  {biomass.id}:")
        print(f"    Functional: {validation['functional']}")
        if validation['functional']:
            print(f"    Objective value: {validation.get('objective_value', 'N/A'):.3f}")
        
        if 'precursor_coverage' in validation:
            print(f"    Precursor coverage:")
            for category, info in validation['precursor_coverage'].items():
                print(f"      {category}: {info['found']}/{info['total']} ({info['percent']:.1f}%)")
    
    # Save model and generate report
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if not args.dry_run and changes_made > 0:
        try:
            write_sbml_model(model, str(output_path))
            print(f"\\n💾 Optimized model saved to: {output_path}")
        except Exception as e:
            print(f"❌ Error saving model: {e}")
            return 1
    
    # Generate report
    report = generate_biomass_report(model, output_path)
    
    print(f"\\n✅ Biomass and NGAM optimization completed!")
    print(f"Changes made: {changes_made}")
    print(f"NGAM present: {report['ngam_present']}")
    
    return 0


if __name__ == "__main__":
    exit(main())