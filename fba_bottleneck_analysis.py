"""
fba_bottleneck_analysis.py

Flux Balance Analysis (FBA) on the E. coli core metabolic model, used to
contextualize how a "bottleneck enzyme" knockout/knockdown affects overall
strain growth and downstream metabolite yield.

This mirrors the kind of "light FBA / pathway-context validation" step
described in Fermelanta's Bioinformatics Internship JD: given a candidate
enzyme identified via docking/screening, use FBA to check whether
engineering it (or a step around it) actually moves the needle on
whole-cell productivity, before investing in wet-lab work.

Model: E. coli "core" metabolic model (Orth, Fleming & Palsson, 2010),
bundled with COBRApy as the standard teaching/benchmark FBA model.

Requirements:
    pip install cobra
"""

import cobra
from cobra.io import load_model
import pandas as pd


def load_ecoli_core():
    """Load the standard E. coli core metabolic model."""
    model = load_model("textbook")
    print(f"Loaded model: {model.id}")
    print(f"  Reactions:   {len(model.reactions)}")
    print(f"  Metabolites: {len(model.metabolites)}")
    print(f"  Genes:       {len(model.genes)}\n")
    return model


def baseline_growth(model):
    """Run FBA under default (aerobic, glucose minimal media) conditions."""
    solution = model.optimize()
    print("=== Baseline FBA (aerobic, glucose minimal media) ===")
    print(f"Growth rate (biomass flux): {solution.objective_value:.4f} h^-1")
    print(f"Status: {solution.status}\n")
    return solution


def top_fluxes(solution, model, n=10):
    """Report the n reactions carrying the largest absolute flux."""
    flux_series = solution.fluxes.abs().sort_values(ascending=False).head(n)
    rows = []
    for rxn_id, flux in flux_series.items():
        rxn = model.reactions.get_by_id(rxn_id)
        rows.append({
            "Reaction": rxn_id,
            "Name": rxn.name,
            "Flux": round(solution.fluxes[rxn_id], 3),
        })
    df = pd.DataFrame(rows)
    print(f"=== Top {n} active reactions by absolute flux ===")
    print(df.to_string(index=False))
    print()
    return df


def knockout_screen(model, reaction_ids):
    """
    Simulate single-reaction knockouts (a stand-in for a rate-limiting /
    'bottleneck' enzyme being non-functional or bottlenecked) and report
    the effect on growth rate relative to baseline.
    """
    baseline = model.optimize().objective_value
    results = []
    for rxn_id in reaction_ids:
        if rxn_id not in model.reactions:
            continue
        with model:  # context manager auto-reverts the knockout after the block
            rxn = model.reactions.get_by_id(rxn_id)
            rxn_name = rxn.name
            model.reactions.get_by_id(rxn_id).knock_out()
            ko_solution = model.optimize()
            ko_growth = ko_solution.objective_value if ko_solution.status == "optimal" else 0.0
        pct_of_baseline = (ko_growth / baseline * 100) if baseline else 0.0
        results.append({
            "Reaction": rxn_id,
            "Name": rxn_name,
            "Baseline growth": round(baseline, 4),
            "Growth after KO": round(ko_growth, 4),
            "% of baseline": round(pct_of_baseline, 1),
        })
    df = pd.DataFrame(results)
    print("=== Bottleneck knockout screen ===")
    print(df.to_string(index=False))
    print()
    return df


def flux_variability(model, reaction_ids):
    """Flux Variability Analysis (FVA) for a shortlist of reactions -
    shows the achievable flux range for each reaction at optimal growth,
    useful for spotting rigid ('bottleneck') vs. flexible steps."""
    from cobra.flux_analysis import flux_variability_analysis
    fva = flux_variability_analysis(model, reaction_list=reaction_ids, fraction_of_optimum=0.9)
    print("=== Flux Variability Analysis (90% of optimal growth) ===")
    print(fva.to_string())
    print()
    return fva


if __name__ == "__main__":
    model = load_ecoli_core()

    baseline_solution = baseline_growth(model)
    top_fluxes(baseline_solution, model, n=10)

    # Candidate "bottleneck" reactions to screen - central carbon metabolism
    # enzymes that are classic rate-limiting steps in E. coli engineering
    # (glycolysis, PPP, and TCA cycle entry points).
    candidate_bottlenecks = [
        "PFK",   # Phosphofructokinase - classic glycolytic bottleneck
        "PYK",   # Pyruvate kinase
        "CS",    # Citrate synthase - TCA cycle entry
        "G6PDH2r",  # Glucose-6-P dehydrogenase - entry to pentose phosphate pathway
        "ENO",   # Enolase
    ]
    ko_df = knockout_screen(model, candidate_bottlenecks)
    fva_df = flux_variability(model, candidate_bottlenecks)

    ko_df.to_csv("knockout_screen_results.csv", index=False)
    fva_df.to_csv("fva_results.csv")
    print("Saved: knockout_screen_results.csv, fva_results.csv")
