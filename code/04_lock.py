#!/usr/bin/env python3
"""Freeze the benchmark and the gate. Run once, commit the output, don't edit."""
import os, hashlib, datetime, textwrap
import pandas as pd

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
ROOT = os.path.expanduser("~/hTYR_Platform")
path = f"{CUR}/benchmark_human_enzyme.csv"

sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
g = pd.read_csv(path)
today = datetime.date.today().strftime("%d-%m-%Y")

n_active = int((g.pAct_median >= 5.5).sum())
lock = f"""# BENCHMARK LOCKED — human tyrosinase

Locked: {today}
File:   data/curated/benchmark_human_enzyme.csv
SHA256: {sha}
N:      {len(g)} compounds, pAct {g.pAct_median.min():.2f}–{g.pAct_median.max():.2f} ({g.pAct_median.max()-g.pAct_median.min():.2f} log span)
Actives (pAct >= 5.5): {n_active}

## Inclusion rules (frozen)
- ChEMBL target CHEMBL1973 (Homo sapiens tyrosinase, UniProt P14679)
- assay_description matched human-enzyme patterns; `unclear` bucket EXCLUDED
- standard_relation '=', units nM, value > 0
- standard_type in IC50/Ki/EC50/Kd
- per compound: median pAct across measurements

## Explicit exclusions
- `unclear` assay descriptions (155 compounds) — provenance unverifiable
- mushroom/fungal and murine-cell assays
- human CELL assays — reserved for validating the full multi-objective stack,
  not the docking layer

## Gate — all three must hold
1. Spearman rho >= 0.30 between predicted score and pAct_median (n={len(g)})
2. Enrichment factor >= 2.0 at top 20%, actives = pAct >= 5.5
3. Rucinol (CHEMBL450195) ranks above kojic acid (CHEMBL287556),
   arbutin (CHEMBL232202) and hydroquinone (CHEMBL537)

Reported but NOT gating: ROC AUC, EF@10%, thiamidol rank, cysteamine rank.
Cysteamine is excluded from gating: its mechanism is not clean active-site
inhibition, so a docking score should not be expected to rank it.

## Pending hand-curation
- Thiamidol (PubChem CID 71543007) — not in ChEMBL; IC50 from Mann et al. 2018
  SMILES: CC(C)C(=O)NC1=NC(=CS1)C2=C(C=C(C=C2)O)O
"""
open(f"{ROOT}/BENCHMARK_LOCKED.md", "w").write(textwrap.dedent(lock))
print(lock)
print(f"wrote {ROOT}/BENCHMARK_LOCKED.md")
