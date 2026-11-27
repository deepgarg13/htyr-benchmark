# Corrections

Errors made during this work, with cause and fix. Published because a benchmark's value depends on
knowing how its numbers were arrived at, including the wrong ones. Two of these produced confident
wrong conclusions that would have redirected the project.

## 1. A calibration objective reported as receptor validation

Until 27-11-2026 this project described the human receptor as validated by a 1.78 A redock of the
co-crystal ligand 0TR. That figure fails as a validation claim in two independent ways.

It is **circular**: 1.78 A is the objective function of the metal-term calibration in
`23_rmsd_select.py`. Scoring parameters were scanned and the set minimising crystal-pose RMSD was
selected (2.09 A baseline, 1.78 A best) -- the log reads `Selected on crystal RMSD alone`.
Reporting the minimised value of one's own selection criterion as evidence of model quality is the
error the pre-registered gate exists to prevent, committed outside the gate's scope.

It is also **in the wrong structure**: that RMSD was measured in the fungal donor crystal 2Y9X, not
in the human model. `19_redock_fixed.py` docks into both receptors; only the 2Y9X block produces an
RMSD. The human block docks hydroquinone, for which no crystallographic pose exists, and reports
scores only.

*Consequence:* **the human receptor has no pose-level validation**, and none is obtainable -- no
crystallographic complex of human tyrosinase with any ligand exists.

*Related tension, now stated in the manuscript:* the metal term used in the best docking
configuration was calibrated in the fungal crystal and applied to the human model, in work arguing
that fungal data do not transfer to human.

*Fix:* the claim is withdrawn from the manuscript, RESULTS.md, the README and the lock file. The
2.09 A baseline and the 1.78 A calibrated value are reported as what they are. The positive control
that would resolve the underlying question -- the identical pipeline against 2Y9X with mushroom
potency data -- is recorded as the most important outstanding experiment.

*Found by tracing the figure back to its log instead of carrying it forward from notes.*

## 2. The frozen gate did not implement the frozen specification

`BENCHMARK_LOCKED.md`, amendment 24-11-2026, replaces condition 1 with condition 1': size-corrected
Spearman >= 0.30. The stated reason is that a null model ranking by heavy-atom count alone scores
rho 0.431 and therefore passes the original condition.

`code/gate.py` continued to compute raw Spearman and test `rho >= 0.30`. The strengthened condition
was written into the specification and never into the code. Attempts 1a–1c were recorded against a
criterion the project had already retired.

A second defect in the same file: condition 4 was skipped silently whenever the decoy was absent
from the scored set, leaving `decoy_ok = True`. An untested condition contributed a pass to the
verdict. A third: rows with missing pAct — the unscored decoy — propagated NaN into the Spearman
calculation, so attempts 2, 2b and 4 recorded `rho = nan` and were failed on the other conditions
alone.

*Fix:* `code/gate_v2.py` implements 1' as specified, drops non-finite rows before correlating,
flags low-coverage subsets explicitly, and distinguishes "condition untested" from "condition
passed". `code/gate.py` is **retained unmodified** for audit. All attempts are re-evaluated under
both gates in `ATTEMPTS.md`.

*Consequence, and it strengthens the conclusion:* under 1', Vina's rho falls from +0.161 to −0.001
and AD4's from +0.086 to −0.243. The apparent docking signal was molecular size throughout.

*Found before publication, by comparing the lock file against the code line by line rather than
trusting that they agreed.*

## 3. Silent subprocess failure reported as a scientific finding (twice)

**`20_metal_visibility.py`** — every docking cell returned an error. The script did not check
return codes, so an empty result matrix was interpreted and reported as "the metal centre is
invisible to the scoring function."

**`23_rmsd_select.py`** — the base scoring terms were built by string-splitting an existing scoring
file, which left a dangling weight and caused every generated scoring file to be rejected. With no
return-code check, the script re-read a stale output file from a previous run, producing 13
identical 2.09 A rows. That near-constant result almost triggered the project's stopping rule.

*Fix:* literal scoring terms instead of string manipulation; `os.remove(out)` before each run;
`raise SystemExit` on any non-zero return code.

*Note on the metal-visibility conclusion:* it was right, but for no reason available at the time.
`19_redock_fixed.py` later returned valid scores -- hydroquinone at -5.3 kcal/mol in the human
receptor with copper, with copper as zinc, and with no metal; 0TR at -6.20 in 2Y9X with and
without copper, at an identical pose. A correct conclusion drawn from a failed computation is
still a failed computation, and the manuscript reports this finding on the strength of the later
runs only.

*Rule adopted:* no script may interpret output it did not confirm was written by a successful
process. Entry 1 above is the same failure in a different costume — a condition that was never
evaluated being counted as satisfied.

## 4. Murine assays classified as human

B16 and other murine cell lines matched the HUMAN regex, producing a false HUMAN_ENZYME
classification in the provenance triage.

*Fix:* separate MURINE class; human regex narrowed and tested against known-murine strings.

## 5. Receptor residue numbering silently broken

pdbfixer renumbered the trimmed receptor from 1 in `10_prep_receptor.py`, and this was not
verified. 15 of 17 key active-site residues were consequently wrong, making the residue-level
interaction features in `38_interaction_fp.py` meaningless.

Caught only because a numbering assertion had been added later.

*Fix:* UniProt numbering recovered by aligning the receptor CA sequence against the P14679 FASTA
(offset +18, 458/458 identities), in `39_fix_numbering.py`. Numbering is now asserted, not assumed.

## 6. Invalid variance comparison

One mushroom subsample was scored 20 times and its cross-validation-seed standard deviation was
presented as comparable to between-subsample variance. These measure different things.

*Fix:* 150 independent subsample draws (`37_percentile.py`). The consistency multiplier in
`RESULTS.md` section 7 comes from the corrected analysis.

## 7. Erroneous censoring note in the lock file

The lock file stated that hydroquinone and arbutin potencies were censored bounds flattened to
exact values. They are exact measured values: Roulier et al., *J Med Chem* 2020, Table 1 reports
hsTYR IC50 of 4,400 uM and 6,500 uM, matching the recorded pAct to within rounding.

*Fix:* struck by dated CORRECTION in `BENCHMARK_LOCKED.md`. The original claim is left visible
rather than deleted.

## 8. Over-strict literature classifier

The paper classifier treated any mention of mushroom tyrosinase as disqualifying, which hid the
most valuable sources: papers assaying both enzymes on the same compound series.

*Fix:* dual-assay papers promoted to highest priority. These supply the paired measurements used in
the transfer analysis.
