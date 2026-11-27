# htyr-benchmark

A locked benchmark and reproducibility harness for **human tyrosinase (hsTYR, UniProt P14679)**
inhibitor ranking, plus five recorded attempts to beat it — four of which failed.

Maintained by Deepanshu Garg. Contributions of human-enzyme data are welcome; see
[Contributing](#contributing).

## Why this exists

Tyrosinase inhibitor discovery has been assayed almost entirely against **mushroom tyrosinase**
(*Agaricus bisporus*, abTYR), a fungal enzyme used as a proxy for the human one. The two share
the binuclear type-3 copper centre and little else: sequence identity in the catalytic domain is
low, the entrance to the active site differs, and compounds that inhibit abTYR potently can be
inert against hsTYR.

The practical consequence is that the field's reference inhibitors rank near the **bottom** of
human-enzyme potency. On the 74-compound benchmark here, arbutin, hydroquinone and kojic acid
place 74th, 73rd and 70th of 74.

This repository does three things:

1. Publishes a **locked, hash-verified benchmark** of human tyrosinase inhibition data, with
   inclusion rules fixed in advance.
2. Publishes a **frozen pass/fail gate** written before any model was scored, so that a method's
   performance cannot be re-specified after the fact.
3. Records **every attempt against that gate, including the failures**, with the corrections log
   for mistakes found during the work.

## The benchmark

| | |
|---|---|
| Compounds | 74 (locked) / 78 evaluated, incl. 3 anchors and 1 decoy |
| Potency span | 4.86 log units |
| Source | ChEMBL, triaged to human-enzyme assays with retained provenance |
| Integrity | SHA256 `5b8d5df2b5b88c6c56cc32060ffaff9e6fb205256235e4b0375f8c0aa5392a6c` |
| Split | dev / holdout, seed 20261123 |

Inclusion rules, amendments and corrections are in [`BENCHMARK_LOCKED.md`](BENCHMARK_LOCKED.md).
Amendments are dated and append-only. Nothing is edited in place.

## The gate

Four conditions, all of which must hold. Frozen in [`code/gate.py`](code/gate.py) before any
method was evaluated.

1'. Size-corrected Spearman rho >= 0.30 — rho of the residual after regressing prediction on
   heavy-atom count, so that a model cannot pass by learning molecular size. Implemented in
   `code/gate_v2.py`; the original `code/gate.py` tested raw rho and is retained unmodified for
   audit (see `CORRECTIONS.md` entry 1)
2. Enrichment factor at 20% >= 2.0
3. Rucinol ranked above kojic acid, arbutin and hydroquinone
4. A known-inactive decoy outside the top 25%

Condition 1 was strengthened from plain rho to size-corrected rho by dated amendment, before
any model passed it.

## Results

| # | Method | Size-corr. rho | EF@20% | Gate |
|---|---|---|---|---|
| 1a | smina, Vina | **−0.001** (raw +0.161) | 1.32 | FAIL |
| 1b | smina, AD4 q=0 | **−0.243** (raw +0.086) | 1.10 | FAIL |
| 1c | smina, AD4 q=+2 | **−0.218** (raw +0.101) | 1.10 | FAIL |
| 2 | best pose, dev | +0.182 | 0.88 | FAIL |
| 4 | calibrated metal term + restricted box, dev | **+0.266** (p=0.061) | 0.59 | FAIL |
| 5 | Random forest, 12 RDKit descriptors, sealed holdout | **+0.467** (perm. p=0.016) | 0.90 | **FAIL** |

**Docking carries no positive ranking information on this target.** Vina's raw correlation of
+0.161 becomes **−0.001** once score is regressed on heavy-atom count; AD4 is significantly
*anti*-correlated at −0.243 (p=0.034). The apparent signal was molecular size. A calibrated metal
term fixed the geometry — closest ligand oxygen to copper moved from 5.74 A to 3.07 A, validated
against crystal-pose RMSD — and did not fix the ranking.

**Ligand-based ML works within scaffold series and collapses across them.** Attempt 5 passed
condition 1 on a random split and failed enrichment. Under a scaffold split (GroupKFold on
Bemis-Murcko frameworks) the same model drops from rho 0.464 to **0.191** — 41% retained.
Structure-based features degrade identically.

That collapse is the central result. A random split made the model look deployable; the scaffold
split shows it cannot generalise to new chemistry. **Generative design against this target is
blocked, not merely suboptimal.**

**The literature has a ceiling.** Extrapolating the scaffold-split learning curve, the entire
1,622-compound aggregated literature reaches only rho 0.347. Interpolating the random-split curve, roughly
**235** compounds measured consistently in one lab would reach rho 0.40. The scaffold-split curve
never reaches it.

**Transfer from mushroom data does not work.** Stacking-based transfer failed. On 17 compounds
with paired human and mushroom measurements, rho is +0.26 (not significant) — the honest reading
is that mushroom potency is *uninformative* about human potency, not that it is anti-correlated.

**The literature gives consistency or breadth, never both.** The human-enzyme subset sits at the
97th percentile of equally sized mushroom subsamples for internal consistency (150 draws), a
multiplier of 8.2–8.5x, at 4.4 compounds per paper against 1.3.

Full numbers and provenance: [`RESULTS.md`](RESULTS.md).

## Receptor construction

No experimental holo structure of full-length human tyrosinase is available. The receptor is
built and validated as follows:

- AlphaFold model of P14679, catalytic domain pLDDT 98.3–98.9
- Binuclear Cu(II) centre transplanted from **PDB 2Y9X** (*A. bisporus* tyrosinase with tropolone) by
  Kabsch superposition of the six-histidine cage, brute-forced over all 36 histidine permutations —
  cage RMSD **0.41 A**
- Trimmed to residues 19–476, protonated with pdbfixer at pH 6.8
- Validation: **none at pose level for the human model** -- no crystallographic complex of human
  tyrosinase exists to redock against. In the donor crystal 2Y9X, 0TR redocks to centroid 1.10 A
  and RMSD 2.09 A, improving to 1.78 A after metal-term calibration; that 1.78 A is the
  calibration objective, not independent validation (`CORRECTIONS.md` entry 1)
- The scoring function does not see the copper: identical scores with the metal present, retyped
  as zinc, and deleted, in both receptors

Note the asymmetry this creates, and it sharpens rather than weakens the thesis: the catalytic
copper cage is conserved closely enough between mushroom and human to superpose at 0.41 A, and
human potencies still do not predict from mushroom potencies. The species difference lives in the
pocket periphery, not the metal site.
- UniProt residue numbering recovered by sequence alignment (offset +18, 458/458 identities)
  after pdbfixer renumbering — see [`code/39_fix_numbering.py`](code/39_fix_numbering.py)

## Reproducing

```bash
conda env create -f environment.yml
conda activate htyr
python code/04_lock.py      # verifies the benchmark hash
python code/41_scaffold_split.py   # the central result
python code/42_scaffold_curve.py
```

Use `code/gate_v2.py` for any new evaluation:

```python
from gate_v2 import run_gate_v2
run_gate_v2(bench, scores, score_col="score", higher_is_better=False)
```

Do not modify either gate. `gate.py` exists only so that the originally recorded numbers can be
reproduced. Record a new attempt in `ATTEMPTS.md` instead.

## Corrections

Eight substantive errors were made during this work, two of which produced confident wrong
conclusions that would have redirected the project, and one of which was a gate that did not
implement its own specification. They are documented with their causes and
fixes in [`CORRECTIONS.md`](CORRECTIONS.md).

That file is part of the record, not an appendix to it. A benchmark whose author does not publish
their own mistakes is not a benchmark.

## What is not claimed

- No compound here is proposed as a development candidate
- No assay was run by the author; all potency values are literature-derived
- The receptor is a validated model, not an experimental structure
- No method in this repository passes the gate

## Contributing

The bottleneck is human-enzyme data, not methods. If you hold human tyrosinase potency
measurements — published or otherwise — adding them to the benchmark raises the ceiling for
everyone working on this target. Open an issue.

Additions follow the inclusion rules in `BENCHMARK_LOCKED.md`, are recorded as dated amendments,
and re-hash the benchmark.

Contact: Deepanshu Garg — deepgarg13599@gmail.com — https://github.com/deepgarg13

## Citation

See [`CITATION.cff`](CITATION.cff).

## Licence

Code: MIT. Benchmark data and documentation: CC-BY-4.0. See [`LICENSE`](LICENSE).
