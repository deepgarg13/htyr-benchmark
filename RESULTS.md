# Results

Every figure here is reproducible from a file in this repository. Attempt-level detail is in
`ATTEMPTS.md`; the benchmark specification and its amendments are in `BENCHMARK_LOCKED.md`.

## 1. The benchmark

| | |
|---|---|
| Locked file | `data/curated/benchmark_human_enzyme.csv` |
| SHA256 | `5b8d5df2b5b88c6c56cc32060ffaff9e6fb205256235e4b0375f8c0aa5392a6c` |
| N | 74 compounds |
| Potency range | pAct 2.19 – 7.05 (4.86 log units) |
| Actives (pAct >= 5.5) | 22 |
| Source target | ChEMBL CHEMBL1973 (*Homo sapiens* tyrosinase, UniProt P14679) |

**Evaluation set.** `data/curated/benchmark_v2.csv`, SHA256 `b1881d29…c176c3de` — 77 scored
compounds plus one unscored decoy. Adds three hand-curated anchors from a single consistent
assay (Mann 2018 Table 1): thiamidol 1.1 uM, 4-hexylresorcinol 94 uM, 4-phenylethylresorcinol
131 uM. The decoy is dimethoxytolylpropyl resorcinol: 0.246 uM against mushroom tyrosinase, **no
inhibition** against the human enzyme. A model that ranks it highly is reproducing mushroom-like
recognition, which is what condition 4 tests.

**Split.** Random, seed 20261123, declared before membership was inspected. dev = 51 (including
the unscored decoy), holdout = 27, SHA256 `e1d9cf8a…ff2bb11e`. The holdout was opened once, for
attempt 5.

**Exclusions.** 155 compounds whose assay descriptions could not be resolved to a species and
preparation; all mushroom, fungal and murine-cell assays; and human *cell* assays, reserved for
validating a full multi-objective stack rather than a docking layer. Cysteamine is present but
excluded from gating — its mechanism is not clean active-site inhibition, so a docking score
should not be expected to rank it.

## 2. The reference inhibitors are at the bottom

Kojic acid, arbutin and hydroquinone — the compounds the field is built on — sit near the floor
of the human-enzyme potency distribution. Arbutin and hydroquinone are recorded at pAct 2.19 and
2.36, matching the exact hsTYR IC50 values of 6,500 uM and 4,400 uM in Roulier et al.,
*J Med Chem* 2020, Table 1.

## 3. Docking carries no positive ranking information on this target

Five configurations. None passed. Under condition 1', the apparent signal is entirely molecular
size:

| configuration | raw rho | size-corrected rho |
|---|---|---|
| Vina | +0.161 (p=0.16) | **−0.001** (p=1.0) |
| AD4, q=0 | +0.086 (p=0.46) | **−0.243** (p=0.034) |
| AD4, q=+2 | +0.101 (p=0.39) | **−0.218** (p=0.059) |
| best pose, dev | −0.237 | +0.182 (p=0.21) |
| calibrated metal + restricted box, dev | −0.264 | **+0.266** (p=0.061) |

Vina's raw correlation of +0.161 becomes −0.001 once score is regressed on heavy-atom count.
AD4 is significantly anti-correlated. Enrichment never exceeded 1.32 and fell below 1.0 on dev.

This is not a weak result but an absent one, and the gate was designed to detect exactly that.
`16_confound_and_chem.py` measured the confound directly: on dev, pAct tracks heavy-atom count at
rho +0.431 (p=1.8e-3, n=50), and a null model ranking by size alone scores EF@20% 1.76. That null
passes the original condition 1 — which is why condition 1 was replaced before any model was
judged against it.

**The physics was fixed; the ranking was not.** A metal term calibrated against crystal-pose RMSD
(`23_rmsd_select.py`, frozen as `models/htyr_frozen.score`) moved the closest ligand oxygen from
5.74 A to 3.07 A of a copper. The best docking attempt that followed reached rho +0.266 — still
under threshold. Making the ligand coordinate the metal correctly did not make the scores rank.

Narrowing to chemically specific poses made matters worse, not better. Requiring a ligand oxygen
within 2.8 A of a copper left only 18 of 50 dev ligands with a qualifying pose, dropped rho to
+0.095, and put the mushroom-potent decoy at rank 1 of 51 — the only attempt to fail condition 4.

## 4. Ligand-based ML works, and fails the gate

Attempt 5: random forest, 500 trees, `min_samples_leaf=2`, on 12 RDKit physicochemical
descriptors. Trained on dev, evaluated once on the sealed holdout.

- size-corrected rho **+0.467**, permutation p = 0.016 over 10,000 label shuffles
- EF@20% **0.90**
- **Gate: FAIL** on condition 2

The only attempt to clear condition 1', and it ranks the top of the list no better than chance.
A model can order a series correctly and still be useless for picking which compounds to make.

## 5. The collapse under scaffold split

This is the central result.

`data/curated/learning_curve.csv` (random split) against `data/curated/scaffold_curve.csv`
(Bemis–Murcko scaffold split), both size-corrected:

| n | random-split rho | scaffold-split rho | scaffolds |
|---|---|---|---|
| ~50 | 0.076 (n=40) | −0.018 | 42 |
| 100 | 0.156 | 0.032 | 73 |
| 200 | — | 0.082 | 122 |
| 400 | 0.452 | 0.206 | 188 |
| 700 | 0.546 (n=650) | 0.245 | 253 |
| 1100 | 0.603 (n=1000) | 0.322 | 316 |
| ~1500 | **0.648** | **0.347** (n=1600) | 385 |

At comparable size, random split reaches rho 0.648 where scaffold split reaches 0.347 — 54%
retained. The same collapse on the single holdout model runs 0.464 to 0.191, 41% retained
(`41_scaffold_split.py`).

**Consequences.**

- The holdout number in section 4 is generous. It came from a random split, which lets structurally
  near-identical compounds fall on both sides.
- Interpolating the random-split curve, size-corrected rho 0.40 arrives at roughly **n = 235**
  consistently measured compounds (0.308 at n=150, 0.417 at n=250).
- The scaffold-split curve **never reaches 0.40**. The entire 1,600-compound aggregated literature,
  spanning 385 scaffolds, tops out at 0.347.
- Generative design against this target is therefore blocked rather than suboptimal: a model that
  cannot rank unseen scaffolds cannot score molecules it has just invented.

Size correction is itself a small-data phenomenon: it moves rho from 0.112 to 0.076 at n=40 and
from 0.6041 to 0.6028 at n=1000. The confound bites hardest in exactly the regime this target
occupies.

## 6. Mushroom data does not transfer

Stacking-based transfer from mushroom to human failed (`34_transfer.py`). On the 17 compounds with
paired measurements, rho is +0.26, not significant. The defensible reading is that mushroom potency
is **uninformative** about human potency — not anti-correlated, which a small positive point
estimate cannot support.

Note the structural asymmetry. The binuclear copper cage used in the human receptor was
transplanted *from mushroom tyrosinase* (PDB 2Y9X) at 0.41 A RMSD. The catalytic centre is
conserved closely enough to superpose almost exactly, and potencies still do not transfer. The
species difference lives in the pocket periphery, not the metal site.

## 7. Consistency and breadth are not simultaneously available

`37_percentile.py`, 150 independent subsample draws: the human-enzyme subset sits at the 97th
percentile of equally sized mushroom subsamples for internal consistency, a multiplier of
8.2–8.5x, at 4.4 compounds per paper against 1.3.

Independently, from the literature rather than from this curation: Roulier et al. 2020 Figure 7 is
the field's complete set of hsTYR inhibitors measured against the isolated enzyme — **25
compounds**. With the Table 1 references and Roulier et al., *Eur J Med Chem* 2022, the published
purified-enzyme SAR totals of order **35–40 compounds**.

Thirty-five consistent compounds, or sixteen hundred inconsistent ones. Neither is enough.

## 8. Receptor construction and validation

| step | result |
|---|---|
| AlphaFold model of P14679, catalytic domain | pLDDT 98.3–98.9 |
| Cu(II) cage transplanted from PDB 2Y9X, Kabsch over all 36 His permutations | RMSD **0.41 A** |
| Trim to residues 19–476, pdbfixer at pH 6.8 | `models/hTYR_receptor.pdb` |
| Redock of co-crystal ligand 0TR (tropolone) | **1.78 A** |
| UniProt numbering recovered by sequence alignment after pdbfixer renumbering | offset +18, 458/458 identities |

The numbering step exists because pdbfixer silently renumbered from 1 and 15 of 17 active-site
residues were consequently wrong. See `CORRECTIONS.md` entry 4.
