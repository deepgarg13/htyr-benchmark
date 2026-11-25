# Gate attempts

Every attempt against the frozen gate, pass or fail, logged before its result was known to be
good. Conditions are defined in `BENCHMARK_LOCKED.md`.

Condition 1' (size-corrected Spearman) is the gating criterion per the 24-11-2026 amendment.
Raw Spearman is reported but does not gate. Attempts 1a–1c were originally recorded against
raw rho only, because `code/gate.py` did not implement 1'; they are re-evaluated here under
`code/gate_v2.py`. See `CORRECTIONS.md` entry 1.

Evaluation set: `data/curated/benchmark_v2.csv` (77 scored + 1 unscored decoy).
Attempts 2–4 were evaluated on dev (51 incl. decoy). Attempt 5 opened the sealed holdout (27).

| # | date | approach | set | n | **1' size-corr rho** | raw rho | 2. EF@20% | 3. anchors | 4. decoy | verdict |
|---|------|----------|-----|---|---------------------|---------|-----------|-----------|----------|---------|
| 1a | 23-11-2026 | smina Vina, neutral ligands, 24 A box, deoxy Cu–Cu 4.36 A | full | 76/78 | **−0.001** (p=1.0) | +0.161 (p=0.16) | 1.32 FAIL | FAIL | n/a | **FAIL** |
| 1b | 23-11-2026 | smina AD4, q=0 | full | 76/78 | **−0.243** (p=0.034) | +0.086 (p=0.46) | 1.10 FAIL | FAIL | n/a | **FAIL** |
| 1c | 23-11-2026 | smina AD4, q=+2 on Cu | full | 76/78 | **−0.218** (p=0.059) | +0.101 (p=0.39) | 1.10 FAIL | FAIL | n/a | **FAIL** |
| 2 | 23-11-2026 | best-scoring pose per ligand | dev | 50 | **+0.182** (p=0.21) | −0.237 (p=0.097) | 0.88 FAIL | FAIL | 14/51 PASS | **FAIL** |
| 2-Cu | 23-11-2026 | best Cu-proximal pose | dev | 30/50 (60%) | +0.238 (p=0.20) † | −0.137 (p=0.47) | 0.88 FAIL | FAIL | 19/51 PASS | **FAIL** |
| 2b | 23-11-2026 | coordinating-pose filter: ligand O within 2.8 A of a Cu | dev | 18/50 (36%) | +0.095 (p=0.71) † | −0.200 (p=0.43) | 0.59 FAIL | FAIL | 1/51 FAIL | **FAIL** |
| 3 | — | no committed score file; not reconstructable | — | — | — | — | — | — | — | — |
| 4 | 24-11-2026 | calibrated metal term (`models/htyr_frozen.score`) + pocket-restricted box | dev | 50 | **+0.266** (p=0.061) | −0.264 (p=0.063) | 0.59 FAIL | FAIL | 18/51 PASS | **FAIL** |
| 5 | 24-11-2026 | RF, 12 RDKit descriptors, dev-trained, SINGLE sealed-holdout evaluation | holdout | 27 | **+0.467** (perm. p=0.016) | — | 0.90 FAIL | — | — | **FAIL** |

† Low coverage. rho is computed over a self-selected subset of ligands that produced a qualifying
pose, and is not comparable to a full-set rho. Flagged automatically by `gate_v2.py`.

Attempt 5 is the only attempt to satisfy condition 1'. It failed condition 2, so the gate did not
open and no generative design was attempted.

Anchor-order detail for attempt 1a — `{rucinol: 58, kojic: 58, arbutin: 33, hydroquinone: 75}`.
Arbutin, which is nearly inactive against the human enzyme, is ranked 25 places better than
rucinol. The species mismatch reproduces inside the gate.
