#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import cobra
import pandas as pd
import os
import json
import numpy as np
from datetime import datetime
import logging
from collections import Counter

# ---------- Utilidades de métricas ----------

def extract_model_name(model_path):
    return os.path.splitext(os.path.basename(model_path))[0]

def load_medium(path):
    """Espera un JSON tipo {"EX_glc__D_e": -10, "EX_o2_e": -1000, ...} con límites inferiores (uptake negativos)."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    # COBRApy espera un dict {exchange_id: bound} -> lo asigna a lower_bound si <0, a ub si >0 según implementación de model.medium
    # Aquí devolvemos tal cual; COBRApy hace el mapeo internamente.
    return data

def safe_len(x):
    try:
        return len(x)
    except Exception:
        return 0

def dead_ends(model):
    prod = Counter()
    cons = Counter()
    for rxn in model.reactions:
        for met, coeff in rxn.metabolites.items():
            if coeff < 0:
                cons[met] += 1
            elif coeff > 0:
                prod[met] += 1
    only_prod = sum((1 for m in prod if cons[m] == 0))
    only_cons = sum((1 for m in cons if prod[m] == 0))
    return only_prod, only_cons

def stoich_sparsity(model):
    m = len(model.metabolites)
    n = len(model.reactions)
    nnz = 0
    for rxn in model.reactions:
        nnz += len(rxn.metabolites)
    total = m * n if m and n else 0
    sparsity = 1.0 - (nnz / total) if total else 1.0
    return sparsity, nnz, total

def degree_means(model):
    met_deg = [len(met.reactions) for met in model.metabolites]
    rxn_deg = [len(rxn.metabolites) for rxn in model.reactions]
    met_mean = float(np.mean(met_deg)) if met_deg else 0.0
    rxn_mean = float(np.mean(rxn_deg)) if rxn_deg else 0.0
    return met_mean, rxn_mean

def gpr_coverage(model):
    return float(np.mean([bool(r.gene_reaction_rule.strip()) for r in model.reactions])) if model.reactions else 0.0

def ec_coverage(model):
    def has_ec(r):
        ann = getattr(r, "annotation", {}) or {}
        ec = ann.get("ec-code") or ann.get("ec") or []
        if isinstance(ec, str):
            ec = [ec]
        return len(ec) > 0
    return float(np.mean([has_ec(r) for r in model.reactions])) if model.reactions else 0.0

def formula_charge_coverage(model):
    has_formula = float(np.mean([bool(getattr(m, 'formula', None)) and getattr(m, 'formula') not in ("", "None")
                                 for m in model.metabolites])) if model.metabolites else 0.0
    has_charge  = float(np.mean([getattr(m, 'charge', None) not in (None, "")
                                 for m in model.metabolites])) if model.metabolites else 0.0
    return has_formula, has_charge

def balance_proxy(model):
    # “Balance” proxy: todas las mets de la rxn tienen formula y charge definidos
    flags = []
    for r in model.reactions:
        mets = r.metabolites.keys()
        ok = all((getattr(m, 'formula', None) not in (None, "", "None") and
                  getattr(m, 'charge',  None) not in (None, "")) for m in mets)
        flags.append(ok)
    return float(np.mean(flags)) if flags else 0.0

def reversible_fraction(model):
    return float(np.mean([r.lower_bound < 0 for r in model.reactions])) if model.reactions else 0.0

def exchanges_summary(model):
    ex = list(model.exchanges)
    if not ex:
        return 0, 0.0, 0.0
    with_lb = float(np.mean([r.lower_bound < 0 for r in ex]))
    with_ub = float(np.mean([r.upper_bound > 0 for r in ex]))
    return len(ex), with_lb, with_ub

def transport_reaction_count(model):
    cnt = 0
    for r in model.reactions:
        comps = {m.compartment for m in r.metabolites}
        if len(comps) >= 2:
            cnt += 1
    return cnt

def currency_load(model, currency_ids=('h_c','h2o_c','atp_c','adp_c','pi_c','ppi_c','nadh_c','nad_c','nadph_c','nadp_c')):
    currency = set(currency_ids)
    if not model.reactions:
        return 0.0
    def has_any(r):
        ids = {m.id for m in r.metabolites}
        return len(ids & currency) > 0
    return float(np.mean([has_any(r) for r in model.reactions]))

def blocked_fraction_fast(model, eps=1e-12):
    # Heurística baratísima por límites estrictos (no sustituye a FVA)
    if not model.reactions:
        return 0.0
    blocked = 0
    for r in model.reactions:
        if abs(r.lower_bound) <= eps and abs(r.upper_bound) <= eps:
            blocked += 1
    return blocked / len(model.reactions)

def find_blocked_optional(model, use_fast):
    if use_fast:
        # Devuelve solo fracción; para conservar compatibilidad con tu campo "blocked_reactions" devolvemos entero aproximado
        frac = blocked_fraction_fast(model)
        return int(round(frac * max(1, len(model.reactions))))
    # Método riguroso (puede ser costoso)
    try:
        from cobra.flux_analysis import find_blocked_reactions
        blocked = find_blocked_reactions(model)
        return len(blocked)
    except Exception:
        # Fallback si falla
        frac = blocked_fraction_fast(model)
        return int(round(frac * max(1, len(model.reactions))))

# ---------- Análisis principal (manteniendo tu interfaz) ----------

def analyze_model(model, biomass_id=None, medium=None, use_fast_blocked=False):
    # Medio (opcional y no destructivo)
    try:
        if medium:
            model.medium = medium
    except Exception as e:
        logging.warning(f"No se pudo aplicar el medio: {e}")

    # FBA (opcionalmente cambiando biomasa)
    try:
        if biomass_id:
            model.objective = biomass_id
        solution = model.optimize()
        objective_value = solution.objective_value if solution.status == 'optimal' else 'Failed'
        fba_status = solution.status
    except Exception as e:
        objective_value = f'Error: {str(e)[:80]}'
        fba_status = 'Error'

    # Métricas rápidas/ampliadas
    only_prod, only_cons = dead_ends(model)
    sparsity, nnz, total = stoich_sparsity(model)
    met_deg_mean, rxn_deg_mean = degree_means(model)
    gpr_cov = gpr_coverage(model)
    ec_cov = ec_coverage(model)
    f_cov, c_cov = formula_charge_coverage(model)
    bal = balance_proxy(model)
    rev_frac = reversible_fraction(model)
    n_ex, ex_lb, ex_ub = exchanges_summary(model)
    n_trans = transport_reaction_count(model)
    currency_frac = currency_load(model)

    blocked_n = 'N/A' if fba_status == 'Error' else find_blocked_optional(model, use_fast_blocked)

    stats = {
        # Núcleo original
        'genes': len(model.genes),
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites),
        'compartments': len(getattr(model, "compartments", {}) or {}),
        'exchange_reactions': n_ex,
        'transport_reactions': n_trans,
        'blocked_reactions': blocked_n,
        'objective_value': objective_value,
        'fba_status': fba_status,
        'biomass_reactions': sum(1 for r in model.reactions if 'biomass' in r.id.lower()),
        'demand_reactions': sum(1 for r in model.reactions if r.id.startswith('DM_')),
        'sink_reactions': sum(1 for r in model.reactions if r.id.startswith('SK_')),

        # Ampliadas
        'gpr_coverage': gpr_cov,                 # 0..1
        'ec_coverage': ec_cov,                   # 0..1
        'met_formula_cov': f_cov,                # 0..1
        'met_charge_cov': c_cov,                 # 0..1
        'balanced_rxn_proxy': bal,               # 0..1
        'reversible_frac': rev_frac,             # 0..1
        'ex_lb_frac': ex_lb,                     # 0..1
        'ex_ub_frac': ex_ub,                     # 0..1
        'dead_end_only_produced': only_prod,     # count
        'dead_end_only_consumed': only_cons,     # count
        'S_sparsity': sparsity,                  # 0..1
        'S_nnz': nnz,                            # int
        'S_total': total,                        # int
        'met_degree_mean': met_deg_mean,         # float
        'rxn_degree_mean': rxn_deg_mean,         # float
        'currency_rxn_frac': currency_frac       # 0..1
    }

    return stats

def difference(a, b):
    return (b - a) if isinstance(a, (int, float, np.number)) and isinstance(b, (int, float, np.number)) else 'N/A'


# ---------- Main ----------
def main():
    parser = argparse.ArgumentParser(description='Compare two SBML metabolic models and export CSV with stats.')
    parser.add_argument('model1', help='Path to first model (.xml/.sbml or .json)')
    parser.add_argument('model2', help='Path to second model (.xml/.sbml or .json)')
    parser.add_argument('-o', '--output', default='reports/', help='Output directory (default: reports/)')
    parser.add_argument('--biomass-id', default=None, help='Reaction ID to set as objective (optional)')
    parser.add_argument('--medium', default=None, help='Path to JSON file with medium (optional)')
    parser.add_argument('--fast-blocked', action='store_true',
                        help='Use fast heuristic for blocked reactions instead of exact method')
    parser.add_argument('--only-comparison', action='store_true',
                        help='Save only the comparison CSV (skip per-model stats CSVs)')

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    os.makedirs(args.output, exist_ok=True)

    try:
        # Carga modelos
        def load_model(path):
            if path.endswith(('.xml', '.sbml')):
                return cobra.io.read_sbml_model(path)
            else:
                return cobra.io.load_json_model(path)

        logging.info(f"Loading model 1: {args.model1}")
        model1 = load_model(args.model1)
        model1_name = os.path.splitext(os.path.basename(args.model1))[0]

        logging.info(f"Loading model 2: {args.model2}")
        model2 = load_model(args.model2)
        model2_name = os.path.splitext(os.path.basename(args.model2))[0]

        # Medio
        medium = None
        if args.medium:
            with open(args.medium, "r", encoding="utf-8") as f:
                medium = json.load(f)

        logging.info("Analyzing models...")
        stats1 = analyze_model(model1, biomass_id=args.biomass_id, medium=medium, use_fast_blocked=args.fast_blocked)
        stats2 = analyze_model(model2, biomass_id=args.biomass_id, medium=medium, use_fast_blocked=args.fast_blocked)

        # Comparación
        metrics = sorted(set(stats1.keys()) | set(stats2.keys()))
        rows = []
        for metric in metrics:
            v1, v2 = stats1.get(metric, 'N/A'), stats2.get(metric, 'N/A')
            rows.append({
                'Metric': metric.replace('_', ' ').title(),
                f'{model1_name}': v1,
                f'{model2_name}': v2,
                'Difference': difference(v1, v2)
            })
        df = pd.DataFrame(rows)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        cmp_file = os.path.join(args.output, f"model_comparison_{model1_name}_vs_{model2_name}_{timestamp}.csv")
        df.to_csv(cmp_file, index=False)

        logging.info(f"Comparison saved to: {cmp_file}")

        # CSVs individuales solo si no se pide --only-comparison
        if not args.only_comparison:
            df1 = pd.DataFrame([stats1]).T.reset_index()
            df1.columns = ['Metric', model1_name]
            m1_file = os.path.join(args.output, f"stats_{model1_name}_{timestamp}.csv")
            df1.to_csv(m1_file, index=False)

            df2 = pd.DataFrame([stats2]).T.reset_index()
            df2.columns = ['Metric', model2_name]
            m2_file = os.path.join(args.output, f"stats_{model2_name}_{timestamp}.csv")
            df2.to_csv(m2_file, index=False)

            logging.info(f"Stats model 1 saved to: {m1_file}")
            logging.info(f"Stats model 2 saved to: {m2_file}")

        print("\nModel Comparison Summary:")
        print(f"Model 1: {model1_name} ({args.model1})")
        print(f"Model 2: {model2_name} ({args.model2})")
        print(f"Comparison CSV: {cmp_file}")
        if not args.only_comparison:
            print(f"Per-model CSVs also saved.")

        print("\nKey differences:")
        for _, row in df.iterrows():
            if row['Difference'] != 'N/A' and row['Difference'] != 0:
                print(f"  {row['Metric']}: {row['Difference']:+}")

    except Exception as e:
        logging.error(f"Error: {str(e)}")
        return 1

    return 0

if __name__ == "__main__":
    exit(main())