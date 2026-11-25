# BENCHMARK LOCKED — human tyrosinase

Locked: 23-11-2026
File:   data/curated/benchmark_human_enzyme.csv
SHA256: 5b8d5df2b5b88c6c56cc32060ffaff9e6fb205256235e4b0375f8c0aa5392a6c
N:      74 compounds, pAct 2.19–7.05 (4.86 log span)
Actives (pAct >= 5.5): 22

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
1. Spearman rho >= 0.30 between predicted score and pAct_median (n=74)
2. Enrichment factor >= 2.0 at top 20%, actives = pAct >= 5.5
3. Rucinol (CHEMBL450195) ranks above kojic acid (CHEMBL287556),
   arbutin (CHEMBL232202) and hydroquinone (CHEMBL537)

Reported but NOT gating: ROC AUC, EF@10%, thiamidol rank, cysteamine rank.
Cysteamine is excluded from gating: its mechanism is not clean active-site
inhibition, so a docking score should not be expected to rank it.

## Pending hand-curation
- Thiamidol (PubChem CID 71543007) — not in ChEMBL; IC50 from Mann et al. 2018
  SMILES: CC(C)C(=O)NC1=NC(=CS1)C2=C(C=C(C=C2)O)O


## AMENDMENT — 23-11-2026
Hand-curated anchors added, as pre-declared above. Source: Mann 2018 Table 1 (hTyr IC50, single consistent assay).
File:   data/curated/benchmark_v2.csv
SHA256: b1881d290c987c3123fd9fefff78f82f459fb8ebb55e49d7663f94dcc176c3de
N:      77 scored compounds (+1 unscored decoy)

Added: Thiamidol (1.1 uM), 4-Hexylresorcinol (94 uM),
       4-Phenylethylresorcinol (131 uM), all hTyr IC50.
Decoy: Dimethoxytolylpropyl resorcinol — mTyr 0.246 uM, hTyr NO INHIBITION.
       Carries no pAct; used only in gate condition 4.

### Gate condition 4 (added)
The decoy must NOT appear in the top 25% of the predicted ranking. A model
that scores it highly is reproducing mushroom-like recognition.

### Data-quality note
ChEMBL records hydroquinone (CHEMBL537) and arbutin (CHEMBL232202) with
standard_relation '=', but Mann 2018 reports both as >4,000 uM. Their pAct
values (2.36, 2.19) are censored bounds flattened to exact values. Harmless
for rank-order use; do not quote them as measured potencies.


## AMENDMENT — 23-11-2026 (dev/holdout split)
After the first gate attempt FAILED on all three conditions, the scored set was
split to protect against overfitting during model iteration.
Random split, seed 20261123, declared before inspecting membership.
  dev     : 51 compounds (incl. the unscored decoy) — iterate freely
  holdout : 27 compounds — SEALED
  holdout SHA256: e1d9cf8acb397ecda4519ed1dd0ac8f8a9b5709b228be210f484a8a7ff2bb11e

Rule: the holdout is opened ONCE, for the first approach that passes the gate on
dev. Gate criteria are unchanged. Every attempt is logged in ATTEMPTS.md,
pass or fail, before its result is known to be good.


## AMENDMENT — 24-11-2026 (gate STRENGTHENED)
Measured on dev: pAct ~ heavy-atom count rho = +0.431 (p=1.8e-3, n=50).
A null model ranking by heavy-atom count alone scores rho 0.431 / EF@20% 1.76,
i.e. it PASSES original condition 1. Condition 1 was therefore too permissive.

Condition 1 is REPLACED (strengthened, never weakened):
  1'. Size-corrected Spearman >= 0.30 — rho between pAct and the residual of the
      predicted score regressed on heavy-atom count. Raw rho is still reported.
Conditions 2, 3 and 4 are unchanged. The null model fails 2 (1.76 < 2.0).
Rationale recorded before any model was evaluated against the new condition.


## CORRECTION — 24-11-2026
The data-quality note above, stating that hydroquinone (CHEMBL537) and arbutin
(CHEMBL232202) are censored bounds flattened to exact values, is WRONG and is
hereby struck. Roulier et al., J Med Chem 2020, Table 1 reports exact hsTYR
IC50 values of 4,400 uM and 6,500 uM respectively. ChEMBL recorded them
correctly; our pAct values (2.36, 2.19) match to within rounding.

## NOTE — world supply of purified-hsTYR data
Roulier et al. 2020 Figure 7 is the field's complete set of hsTYR inhibitors
measured against the isolated enzyme: 25 compounds. With Table 1 references and
Roulier et al., Eur J Med Chem 2022 (compounds 10 and 23, Ki 0.55 and 0.25 uM),
the published purified-enzyme SAR is of order 35-40 compounds. This independently
confirms the data-scarcity finding reached from ChEMBL and BindingDB.


## AMENDMENT — 25-11-2026 (gate implementation corrected)
The 24-11-2026 amendment replaced condition 1 with condition 1' (size-corrected Spearman).
code/gate.py was never updated and continued to test raw Spearman, i.e. the retired condition.
Condition 4 was also skipped silently when the decoy was absent from the scored set, and missing
pAct values propagated NaN into the correlation.

code/gate_v2.py is added, implementing 1' as specified, dropping non-finite rows before
correlating, flagging low-coverage subsets, and distinguishing an untested condition from a
passed one. code/gate.py is RETAINED UNMODIFIED for audit; it must not be used for new
evaluations. The gate criteria themselves are unchanged — only their implementation is fixed.

All attempts re-evaluated under both gates; see ATTEMPTS.md. Under 1', Vina falls from raw
+0.161 to -0.001 and AD4 q0 from +0.086 to -0.243. Attempt 5 (+0.467 holdout) is unaffected:
25_holdout.py computed 1' inline and never imported gate.py.
