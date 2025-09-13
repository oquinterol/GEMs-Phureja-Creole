#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import cobra


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Smoke FBA: load model, optimize, log objective and basic stats")
    p.add_argument("--model", required=True, help="Path to SBML model")
    p.add_argument("--out", required=True, help="Output log path")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    model = cobra.io.read_sbml_model(args.model)
    # Prefer GLPK if available
    try:
        model.solver = "glpk"
    except Exception:
        pass
    res = model.optimize()
    # Collect stats
    info = {
        "status": str(res.status),
        "objective_value": float(res.objective_value) if res.objective_value is not None else None,
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
        "n_genes": len(model.genes),
        "objective_reaction": getattr(model.objective, "expression", None).__str__() if model.objective is not None else None,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("SMOKE FBA\n")
        fh.write(json.dumps(info, indent=2))
        fh.write("\n")
    print(f"Wrote smoke log to {out}")


if __name__ == "__main__":
    main()

