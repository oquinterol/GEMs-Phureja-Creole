#!/usr/bin/env python3
"""
Mass and charge balance correction script for SBML metabolic models.
Identifies and attempts to fix stoichiometric imbalances.

Author: Claude Code
Date: 2025-09-13
"""

import os
import sys
import json
import logging
import argparse
import cobra
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import re

class ModelBalancer:
    """Main class for balancing metabolic model reactions."""

    def __init__(self, model_path: str, output_path: str):
        self.model_path = Path(model_path)
        self.output_path = Path(output_path)

        # Load model
        self.model = cobra.io.read_sbml_model(str(self.model_path))

        # Setup logging
        self.setup_logging()

        # Balance tracking
        self.balance_report = {
            'total_reactions': len(self.model.reactions),
            'mass_imbalanced': [],
            'charge_imbalanced': [],
            'fixed_mass': [],
            'fixed_charge': [],
            'unfixable': [],
            'summary': {}
        }

        # Common molecular formulas for metabolites
        self.formula_patterns = {
            'cpd00001': 'H2O',      # Water
            'cpd00002': 'C21H28N7O17P3',  # ATP
            'cpd00008': 'C21H27N7O17P3',  # ADP
            'cpd00009': 'H3O4P',     # Phosphate
            'cpd00067': 'H',         # Proton
            'cpd00003': 'C21H27N7O17P3',  # NAD
            'cpd00004': 'C21H26N7O17P3',  # NADH
            'cpd00006': 'C21H27N7O17P3',  # NADP
            'cpd00005': 'C21H28N7O18P3',  # NADPH
            'cpd00011': 'CO2',       # Carbon dioxide
            'cpd00007': 'O2',        # Oxygen
        }

    def setup_logging(self):
        """Setup comprehensive logging."""
        log_file = self.output_path.parent / f"balance_{int(time.time())}.log"

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Starting balance correction for {self.model_path}")

    def parse_formula(self, formula: str) -> Dict[str, int]:
        """Parse molecular formula into element counts."""
        if not formula or formula == "":
            return {}

        elements = defaultdict(int)

        # Remove charges and handle special cases
        formula = formula.replace('+', '').replace('-', '')

        # Pattern to match element symbols and counts
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)

        for element, count in matches:
            count = int(count) if count else 1
            elements[element] += count

        return dict(elements)

    def get_metabolite_formula(self, metabolite) -> Optional[str]:
        """Get metabolite formula from various sources."""
        # Try model annotation first
        if hasattr(metabolite, 'formula') and metabolite.formula:
            return metabolite.formula

        # Try known patterns
        compound_id = metabolite.id.split('_')[0]
        if compound_id in self.formula_patterns:
            return self.formula_patterns[compound_id]

        # Try to extract from name or annotation
        if hasattr(metabolite, 'annotation'):
            for key, value in metabolite.annotation.items():
                if 'formula' in key.lower():
                    return value

        return None

    def calculate_mass_balance(self, reaction) -> Dict[str, float]:
        """Calculate elemental mass balance for a reaction."""
        element_balance = defaultdict(float)

        for metabolite, coefficient in reaction.metabolites.items():
            formula = self.get_metabolite_formula(metabolite)
            if not formula:
                # Can't balance without formula
                return {'unknown': 1.0}

            elements = self.parse_formula(formula)
            for element, count in elements.items():
                element_balance[element] += coefficient * count

        return dict(element_balance)

    def calculate_charge_balance(self, reaction) -> float:
        """Calculate charge balance for a reaction."""
        total_charge = 0.0

        for metabolite, coefficient in reaction.metabolites.items():
            charge = getattr(metabolite, 'charge', 0) or 0
            total_charge += coefficient * charge

        return total_charge

    def is_transport_reaction(self, reaction) -> bool:
        """Check if reaction is a transport reaction."""
        # Transport reactions move metabolites between compartments
        compartments = set()
        for metabolite in reaction.metabolites:
            comp = metabolite.compartment
            compartments.add(comp)

        # If more than one compartment, likely transport
        return len(compartments) > 1

    def is_exchange_reaction(self, reaction) -> bool:
        """Check if reaction is an exchange/boundary reaction."""
        return (reaction.id.startswith('EX_') or
                reaction.id.startswith('DM_') or
                reaction.id.startswith('SK_') or
                len(reaction.metabolites) == 1)

    def fix_proton_balance(self, reaction, element_balance: Dict[str, float]) -> bool:
        """Attempt to fix hydrogen imbalance by adding/removing protons."""
        h_imbalance = element_balance.get('H', 0)

        if abs(h_imbalance) < 0.01:  # Already balanced
            return True

        # Find proton in appropriate compartment
        proton_met = None
        main_compartment = None

        # Determine main compartment
        for metabolite in reaction.metabolites:
            if not metabolite.id.startswith('cpd00067'):  # Not a proton
                main_compartment = metabolite.compartment
                break

        # Find or create proton
        for metabolite in self.model.metabolites:
            if (metabolite.id.startswith('cpd00067') and
                metabolite.compartment == main_compartment):
                proton_met = metabolite
                break

        if proton_met:
            # Add protons to balance
            if proton_met in reaction.metabolites:
                reaction.metabolites[proton_met] -= h_imbalance
            else:
                reaction.add_metabolites({proton_met: -h_imbalance})

            self.logger.info(f"Fixed H balance in {reaction.id} by adjusting protons by {-h_imbalance}")
            return True

        return False

    def fix_water_balance(self, reaction, element_balance: Dict[str, float]) -> bool:
        """Attempt to fix oxygen imbalance by adding/removing water."""
        o_imbalance = element_balance.get('O', 0)
        h_imbalance = element_balance.get('H', 0)

        # Water has H2O, so O imbalance should be accompanied by 2*H imbalance
        if abs(o_imbalance) < 0.01:
            return True

        if abs(h_imbalance - 2 * o_imbalance) > 0.1:
            # Imbalance not consistent with water
            return False

        # Find water in appropriate compartment
        water_met = None
        main_compartment = None

        for metabolite in reaction.metabolites:
            main_compartment = metabolite.compartment
            break

        for metabolite in self.model.metabolites:
            if (metabolite.id.startswith('cpd00001') and
                metabolite.compartment == main_compartment):
                water_met = metabolite
                break

        if water_met:
            if water_met in reaction.metabolites:
                reaction.metabolites[water_met] -= o_imbalance
            else:
                reaction.add_metabolites({water_met: -o_imbalance})

            self.logger.info(f"Fixed O/H balance in {reaction.id} by adjusting water by {-o_imbalance}")
            return True

        return False

    def balance_reaction_mass(self, reaction) -> bool:
        """Attempt to balance a single reaction for mass."""
        # Skip certain reaction types
        if self.is_exchange_reaction(reaction) or self.is_transport_reaction(reaction):
            return True

        element_balance = self.calculate_mass_balance(reaction)

        # Check if already balanced
        max_imbalance = max(abs(val) for val in element_balance.values()) if element_balance else 0
        if max_imbalance < 0.01:
            return True

        # Try common fixes
        original_metabolites = dict(reaction.metabolites)

        # Try to fix hydrogen balance with protons
        if 'H' in element_balance and abs(element_balance['H']) > 0.01:
            if self.fix_proton_balance(reaction, element_balance):
                # Recalculate balance
                element_balance = self.calculate_mass_balance(reaction)

        # Try to fix oxygen balance with water
        if 'O' in element_balance and abs(element_balance['O']) > 0.01:
            if self.fix_water_balance(reaction, element_balance):
                element_balance = self.calculate_mass_balance(reaction)

        # Check if fixed
        max_imbalance = max(abs(val) for val in element_balance.values()) if element_balance else 0
        if max_imbalance < 0.01:
            self.balance_report['fixed_mass'].append(reaction.id)
            return True
        else:
            # Revert changes if not fixed
            reaction.metabolites.clear()
            reaction.add_metabolites(original_metabolites)
            return False

    def balance_reaction_charge(self, reaction) -> bool:
        """Attempt to balance a single reaction for charge."""
        if self.is_exchange_reaction(reaction):
            return True

        charge_balance = self.calculate_charge_balance(reaction)

        if abs(charge_balance) < 0.01:
            return True

        # Try to fix by adding protons
        main_compartment = None
        for metabolite in reaction.metabolites:
            main_compartment = metabolite.compartment
            break

        # Find proton
        proton_met = None
        for metabolite in self.model.metabolites:
            if (metabolite.id.startswith('cpd00067') and
                metabolite.compartment == main_compartment):
                proton_met = metabolite
                break

        if proton_met:
            if proton_met in reaction.metabolites:
                reaction.metabolites[proton_met] -= charge_balance
            else:
                reaction.add_metabolites({proton_met: -charge_balance})

            self.logger.info(f"Fixed charge balance in {reaction.id} by adjusting protons by {-charge_balance}")
            self.balance_report['fixed_charge'].append(reaction.id)
            return True

        return False

    def analyze_balance(self):
        """Analyze current balance status of all reactions."""
        self.logger.info("Analyzing mass and charge balance...")

        for reaction in self.model.reactions:
            # Mass balance
            element_balance = self.calculate_mass_balance(reaction)
            if element_balance and 'unknown' not in element_balance:
                max_imbalance = max(abs(val) for val in element_balance.values())
                if max_imbalance > 0.01:
                    self.balance_report['mass_imbalanced'].append({
                        'id': reaction.id,
                        'imbalance': element_balance
                    })

            # Charge balance
            charge_balance = self.calculate_charge_balance(reaction)
            if abs(charge_balance) > 0.01:
                self.balance_report['charge_imbalanced'].append({
                    'id': reaction.id,
                    'imbalance': charge_balance
                })

    def fix_all_reactions(self):
        """Attempt to fix all imbalanced reactions."""
        self.logger.info("Attempting to fix mass imbalances...")

        for reaction_info in self.balance_report['mass_imbalanced']:
            reaction = self.model.reactions.get_by_id(reaction_info['id'])
            if not self.balance_reaction_mass(reaction):
                self.balance_report['unfixable'].append({
                    'id': reaction.id,
                    'type': 'mass',
                    'imbalance': reaction_info['imbalance']
                })

        self.logger.info("Attempting to fix charge imbalances...")

        for reaction_info in self.balance_report['charge_imbalanced']:
            reaction = self.model.reactions.get_by_id(reaction_info['id'])
            if not self.balance_reaction_charge(reaction):
                self.balance_report['unfixable'].append({
                    'id': reaction.id,
                    'type': 'charge',
                    'imbalance': reaction_info['imbalance']
                })

    def generate_summary(self):
        """Generate balance summary statistics."""
        self.balance_report['summary'] = {
            'total_reactions': self.balance_report['total_reactions'],
            'mass_imbalanced_initial': len(self.balance_report['mass_imbalanced']),
            'charge_imbalanced_initial': len(self.balance_report['charge_imbalanced']),
            'mass_fixed': len(self.balance_report['fixed_mass']),
            'charge_fixed': len(self.balance_report['fixed_charge']),
            'unfixable': len(self.balance_report['unfixable']),
            'mass_balance_rate': (self.balance_report['total_reactions'] -
                                 len(self.balance_report['mass_imbalanced']) +
                                 len(self.balance_report['fixed_mass'])) / self.balance_report['total_reactions'],
            'charge_balance_rate': (self.balance_report['total_reactions'] -
                                   len(self.balance_report['charge_imbalanced']) +
                                   len(self.balance_report['fixed_charge'])) / self.balance_report['total_reactions']
        }

    def save_results(self):
        """Save balanced model and report."""
        # Save model
        cobra.io.write_sbml_model(self.model, str(self.output_path))

        # Save report
        report_path = self.output_path.parent / f"{self.output_path.stem}_balance_report.json"
        with open(report_path, 'w') as f:
            json.dump(self.balance_report, f, indent=2)

        # Save summary
        summary_path = self.output_path.parent / f"{self.output_path.stem}_balance_summary.txt"
        with open(summary_path, 'w') as f:
            f.write("MASS AND CHARGE BALANCE REPORT\n")
            f.write("="*50 + "\n\n")
            f.write(f"Total reactions: {self.balance_report['summary']['total_reactions']}\n")
            f.write(f"Mass imbalanced (initial): {self.balance_report['summary']['mass_imbalanced_initial']}\n")
            f.write(f"Charge imbalanced (initial): {self.balance_report['summary']['charge_imbalanced_initial']}\n")
            f.write(f"Mass balance fixed: {self.balance_report['summary']['mass_fixed']}\n")
            f.write(f"Charge balance fixed: {self.balance_report['summary']['charge_fixed']}\n")
            f.write(f"Unfixable reactions: {self.balance_report['summary']['unfixable']}\n")
            f.write(f"Final mass balance rate: {self.balance_report['summary']['mass_balance_rate']:.1%}\n")
            f.write(f"Final charge balance rate: {self.balance_report['summary']['charge_balance_rate']:.1%}\n")

        self.logger.info(f"Results saved to {self.output_path}")
        self.logger.info(f"Report saved to {report_path}")
        self.logger.info(f"Summary saved to {summary_path}")

    def run_balancing(self):
        """Run complete balancing pipeline."""
        try:
            self.logger.info("Starting mass and charge balancing pipeline")

            # Phase 1: Analyze current state
            self.analyze_balance()

            # Phase 2: Attempt fixes
            self.fix_all_reactions()

            # Phase 3: Generate summary
            self.generate_summary()

            # Phase 4: Save results
            self.save_results()

            # Log results
            summary = self.balance_report['summary']
            self.logger.info("="*60)
            self.logger.info("BALANCING COMPLETED")
            self.logger.info(f"Mass balance rate: {summary['mass_balance_rate']:.1%}")
            self.logger.info(f"Charge balance rate: {summary['charge_balance_rate']:.1%}")
            self.logger.info(f"Reactions fixed: {summary['mass_fixed'] + summary['charge_fixed']}")
            self.logger.info(f"Unfixable reactions: {summary['unfixable']}")
            self.logger.info("="*60)

        except Exception as e:
            self.logger.error(f"Balancing failed: {e}")
            raise


def main():
    parser = argparse.ArgumentParser(description="Mass and charge balance correction")
    parser.add_argument("--model", required=True, help="Input SBML model path")
    parser.add_argument("--output", required=True, help="Output SBML model path")

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.model).exists():
        print(f"Error: Input model {args.model} does not exist")
        sys.exit(1)

    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # Run balancing
    balancer = ModelBalancer(args.model, args.output)
    balancer.run_balancing()


if __name__ == "__main__":
    main()