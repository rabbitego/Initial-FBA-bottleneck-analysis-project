# Flux Balance Analysis: Enzyme Bottleneck Validation in E. coli

This project demonstrates a lightweight Flux Balance Analysis (FBA) workflow for evaluating candidate enzyme bottlenecks in central carbon metabolism. The goal is to connect a biological enzyme target or reaction candidate with whole-cell metabolic consequences, rather than treating the enzyme in isolation.

The workflow mirrors the kind of pathway-context validation expected in a computational biology or enzyme engineering setting: use a metabolic model to test whether a candidate bottleneck actually constrains growth or production realistically.

## Project Objective

Computational enzyme discovery tools such as docking, structure prediction, and scoring can identify enzyme variants with improved binding affinity or predicted activity. However, these tools do not show whether a change in that enzyme will improve the final strain-level phenotype.

This project closes that gap by using FBA to ask a system-level question:

> Given the E. coli core metabolic network and a default aerobic glucose-minimal growth condition, how much does disabling each candidate reaction reduce the achievable growth rate, and how flexible is that part of the network?

## Workflow

1. Load the standard E. coli core metabolic model bundled with COBRApy.
2. Run a baseline FBA optimization to estimate the optimal growth rate.
3. Rank the top active reactions by absolute flux magnitude.
4. Perform a knockout screen over five candidate reactions:
   - PFK
   - PYK
   - CS
   - G6PDH2r
   - ENO
5. Perform Flux Variability Analysis (FVA) at 90% of optimal growth for the same reactions.

## Biological Interpretation

The knockout screen estimates how much growth remains after each reaction is removed. FVA then reveals whether the reaction is tightly constrained or flexible:

- A narrow flux range indicates a rigid, highly constrained reaction.
- A wide flux range indicates a flexible step that the network can compensate around.

This creates a simple decision framework for enzyme target prioritization:

- Essential reactions are poor bottleneck engineering targets because they are required for viability.
- Strong bottlenecks are reactions whose removal causes a measurable but not catastrophic loss in growth.
- Non-limiting reactions are weak engineering targets because the network can reroute around them.

## Example Results

| Enzyme | Growth after knockout | Interpretation |
|---|---:|---|
| Citrate synthase (CS) | 0% of baseline | Essential reaction; no apparent bypass exists |
| Enolase (ENO) | ~0% of baseline | Essential core glycolytic step |
| Phosphofructokinase (PFK) | 80.6% of baseline | Genuine bottleneck with partial compensation |
| Pyruvate kinase (PYK) | 99.0% of baseline | Non-limiting, highly compensable |
| Glucose-6-P dehydrogenase (G6PDH2r) | 98.8% of baseline | Non-limiting, highly compensable |

### Takeaway

Among the tested reactions, PFK behaves most like a meaningful bottleneck target for engineering or pathway optimization. CS and ENO are essential steps, while PYK and G6PDH2r appear non-limiting because the network can reroute around them.

## Outputs

This project writes two analysis tables:

- knockout_screen_results.csv
- fva_results.csv

## Repository Structure

```text
.
├── fba_bottleneck_analysis.py
├── knockout_screen_results.csv
├── fva_results.csv
├── README.md
└── requirements.txt
```

## Setup and Run

```bash
pip install -r requirements.txt
python fba_bottleneck_analysis.py
```

## Scope and Honest Limitations

This project uses the standard E. coli core metabolic model with 95 reactions. It is intended as a teaching and demonstration FBA workflow rather than a full genome-scale production model. A full strain-specific or genome-scale model would provide more realistic engineering conclusions, but the core model is appropriate for an internship-style pathway-context analysis.

