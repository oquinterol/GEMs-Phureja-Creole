#!/usr/bin/env python3
"""
Create unified base model from drafts and curated versions.

Takes the best components from multiple models to create a unified base model:
- Uses the most complete/curated model as base
- Optionally merges reactions/metabolites from drafts if they add value
- Validates the unified model before saving

Usage:
  python scripts/create_unified_model.py \
    --base models/curated/creole_with_gprs_annot.xml \
    --drafts models/drafts/creole_kbase.xml models/drafts/Botara.xml \
    --out models/current/creole_unified.xml
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Set
import cobra
from cobra.io import read_sbml_model, write_sbml_model


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create unified model from drafts and curated versions")
    p.add_argument("--base", required=True, help="Base model path (most complete/curated)")
    p.add_argument("--drafts", nargs="+", help="Draft model paths to merge components from")
    p.add_argument("--out", required=True, help="Output unified model path")
    p.add_argument("--merge-reactions", action="store_true", 
                   help="Merge unique reactions from drafts (default: use base only)")
    p.add_argument("--merge-metabolites", action="store_true",
                   help="Merge unique metabolites from drafts (default: use base only)")
    p.add_argument("--validate", action="store_true", default=True,
                   help="Validate model after unification (default: True)")
    return p.parse_args()


def get_model_stats(model: cobra.Model) -> dict:
    """Get basic statistics from a model."""
    return {
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
        "n_genes": len(model.genes),
        "model_id": model.id,
        "reactions": [r.id for r in model.reactions],
        "metabolites": [m.id for m in model.metabolites],
        "genes": [g.id for g in model.genes]
    }


def merge_unique_reactions(base_model: cobra.Model, draft_models: List[cobra.Model]) -> int:
    """Merge unique reactions from draft models into base model."""
    base_rxn_ids = {r.id for r in base_model.reactions}
    added_count = 0
    
    for draft in draft_models:
        for reaction in draft.reactions:
            if reaction.id not in base_rxn_ids:
                # Check if metabolites exist in base model
                missing_mets = []
                for met in reaction.metabolites:
                    if met.id not in [m.id for m in base_model.metabolites]:
                        missing_mets.append(met)
                
                # Add missing metabolites first
                for met in missing_mets:
                    try:
                        base_model.add_metabolites([met.copy()])
                    except Exception as e:
                        print(f"Warning: Could not add metabolite {met.id}: {e}")
                        continue
                
                # Add the reaction
                try:
                    base_model.add_reactions([reaction.copy()])
                    added_count += 1
                    base_rxn_ids.add(reaction.id)
                    print(f"Added reaction: {reaction.id}")
                except Exception as e:
                    print(f"Warning: Could not add reaction {reaction.id}: {e}")
    
    return added_count


def merge_unique_metabolites(base_model: cobra.Model, draft_models: List[cobra.Model]) -> int:
    """Merge unique metabolites from draft models into base model."""
    base_met_ids = {m.id for m in base_model.metabolites}
    added_count = 0
    
    for draft in draft_models:
        for metabolite in draft.metabolites:
            if metabolite.id not in base_met_ids:
                try:
                    base_model.add_metabolites([metabolite.copy()])
                    added_count += 1
                    base_met_ids.add(metabolite.id)
                    print(f"Added metabolite: {metabolite.id}")
                except Exception as e:
                    print(f"Warning: Could not add metabolite {metabolite.id}: {e}")
    
    return added_count


def validate_model(model: cobra.Model) -> dict:
    """Basic model validation."""
    try:
        # Try FBA
        solution = model.optimize()
        
        # Check for basic consistency
        errors = []
        warnings = []
        
        # Check for orphaned metabolites
        orphaned_mets = []
        for met in model.metabolites:
            if len(met.reactions) == 0:
                orphaned_mets.append(met.id)
        
        if orphaned_mets:
            warnings.append(f"Found {len(orphaned_mets)} orphaned metabolites")
        
        # Check for reactions without GPRs (if genes exist)
        if len(model.genes) > 0:
            no_gpr_rxns = [r.id for r in model.reactions if not r.gene_reaction_rule]
            if no_gpr_rxns:
                warnings.append(f"Found {len(no_gpr_rxns)} reactions without GPRs")
        
        validation_result = {
            "valid": solution.status == "optimal",
            "objective_value": float(solution.objective_value) if solution.objective_value else None,
            "status": str(solution.status),
            "errors": errors,
            "warnings": warnings,
            "orphaned_metabolites": len(orphaned_mets),
            "reactions_without_gpr": len(no_gpr_rxns) if len(model.genes) > 0 else 0
        }
        
        return validation_result
        
    except Exception as e:
        return {
            "valid": False,
            "error": str(e),
            "errors": [str(e)],
            "warnings": []
        }


def main() -> None:
    args = parse_args()
    
    print(f"Creating unified model from base: {args.base}")
    
    # Load base model
    base_model = read_sbml_model(args.base)
    print(f"Base model: {base_model.id}")
    base_stats = get_model_stats(base_model)
    print(f"Base stats: {base_stats['n_reactions']} reactions, {base_stats['n_metabolites']} metabolites, {base_stats['n_genes']} genes")
    
    # Load draft models if provided
    draft_models = []
    if args.drafts:
        for draft_path in args.drafts:
            try:
                draft = read_sbml_model(draft_path)
                draft_models.append(draft)
                draft_stats = get_model_stats(draft)
                print(f"Draft {draft.id}: {draft_stats['n_reactions']} reactions, {draft_stats['n_metabolites']} metabolites, {draft_stats['n_genes']} genes")
            except Exception as e:
                print(f"Warning: Could not load draft model {draft_path}: {e}")
    
    # Create unified model (start with copy of base)
    unified_model = base_model.copy()
    unified_model.id = "CreolePotato_Unified"
    unified_model.name = "Unified Creole Potato Metabolic Model"
    
    # Merge components if requested
    added_reactions = 0
    added_metabolites = 0
    
    if args.merge_reactions and draft_models:
        print("Merging unique reactions from drafts...")
        added_reactions = merge_unique_reactions(unified_model, draft_models)
        print(f"Added {added_reactions} unique reactions")
    
    if args.merge_metabolites and draft_models:
        print("Merging unique metabolites from drafts...")
        added_metabolites = merge_unique_metabolites(unified_model, draft_models)
        print(f"Added {added_metabolites} unique metabolites")
    
    # Final stats
    final_stats = get_model_stats(unified_model)
    print(f"Unified model stats: {final_stats['n_reactions']} reactions, {final_stats['n_metabolites']} metabolites, {final_stats['n_genes']} genes")
    
    # Validate if requested
    if args.validate:
        print("Validating unified model...")
        validation = validate_model(unified_model)
        print(f"Validation result: {'PASS' if validation['valid'] else 'FAIL'}")
        print(f"FBA status: {validation.get('status', 'unknown')}")
        if validation.get('objective_value'):
            print(f"Objective value: {validation['objective_value']:.6f}")
        
        if validation.get('warnings'):
            print("Warnings:")
            for warning in validation['warnings']:
                print(f"  - {warning}")
        
        if validation.get('errors'):
            print("Errors:")
            for error in validation['errors']:
                print(f"  - {error}")
    
    # Save unified model
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        write_sbml_model(unified_model, str(out_path))
        print(f"Unified model saved to: {out_path}")
        
        # Save summary
        summary = {
            "base_model": args.base,
            "draft_models": args.drafts or [],
            "base_stats": base_stats,
            "final_stats": final_stats,
            "added_reactions": added_reactions,
            "added_metabolites": added_metabolites,
            "validation": validation if args.validate else None
        }
        
        summary_path = out_path.with_suffix('.summary.json')
        with summary_path.open('w') as f:
            json.dump(summary, f, indent=2)
        print(f"Summary saved to: {summary_path}")
        
    except Exception as e:
        print(f"Error saving unified model: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())