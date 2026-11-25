"""Gate v2 — implements condition 1' from BENCHMARK_LOCKED.md (amendment 24-11-2026).

code/gate.py is retained unmodified for audit. Differences here:
  1. Condition 1 is SIZE-CORRECTED Spearman, as the lock file specifies.
     gate.py tested raw Spearman, the condition the amendment retired.
     Raw rho is still reported alongside, per the amendment.
  2. Rows with missing pAct_median (the unscored decoy) are dropped before
     computing rho and EF. In gate.py they produced rho = nan.
  3. Condition 4 prints SKIPPED when the decoy is absent from the scored set,
     and marks the verdict INCOMPLETE rather than passing untested.
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ANCHOR = dict(rucinol="CHEMBL450195", kojic="CHEMBL287556",
              arbutin="CHEMBL232202", hydroquinone="CHEMBL537")
DECOY_ROLE = "decoy_mushroom_potent"
ACTIVE_CUT = 5.5


def _heavy_atoms(smiles):
    from rdkit import Chem
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        out.append(np.nan if m is None else m.GetNumHeavyAtoms())
    return np.asarray(out, dtype=float)


def _size_corrected_rho(score, pact, heavy):
    """rho between pAct and the residual of score regressed on heavy-atom count."""
    ok = np.isfinite(score) & np.isfinite(pact) & np.isfinite(heavy)
    if ok.sum() < 5:
        return np.nan, np.nan, int(ok.sum())
    A = np.column_stack([np.ones(ok.sum()), heavy[ok]])
    beta, *_ = np.linalg.lstsq(A, score[ok], rcond=None)
    resid = score[ok] - A @ beta
    rho, p = spearmanr(resid, pact[ok])
    return rho, p, int(ok.sum())


def run_gate_v2(bench, scores, score_col="score", higher_is_better=True,
                decoy_expected=True):
    d = bench.merge(scores, on="molecule_chembl_id", how="inner")
    d = d.assign(_s=d[score_col] if higher_is_better else -d[score_col])

    # decoy is handled separately; it carries no pAct
    is_decoy = (d.role == DECOY_ROLE) if "role" in d else pd.Series(False, index=d.index)
    scored = d[~is_decoy & d.pAct_median.notna()].copy()

    if "heavy" in scored:
        heavy = scored.heavy.to_numpy(dtype=float)
    elif "smiles" in scored:
        heavy = _heavy_atoms(scored.smiles)
    else:
        raise ValueError("need a 'heavy' column or 'smiles' to compute it")

    _fin = np.isfinite(scored._s.to_numpy(dtype=float)) & scored.pAct_median.notna().to_numpy()
    raw_rho, raw_p = (spearmanr(scored._s[_fin], scored.pAct_median[_fin])
                      if _fin.sum() >= 5 else (np.nan, np.nan))
    sc_rho, sc_p, n_sc = _size_corrected_rho(
        scored._s.to_numpy(dtype=float), scored.pAct_median.to_numpy(dtype=float), heavy)

    ranked = scored.sort_values("_s", ascending=False)
    k = max(1, int(round(0.20 * len(ranked))))
    base = (ranked.pAct_median >= ACTIVE_CUT).mean()
    ef = ((ranked.head(k).pAct_median >= ACTIVE_CUT).sum() / k / base) if base else np.nan

    rank = {n: int((ranked._s > ranked.loc[ranked.molecule_chembl_id == c, "_s"].iloc[0]).sum()) + 1
            for n, c in ANCHOR.items() if (ranked.molecule_chembl_id == c).any()}
    ordering = all(rank.get("rucinol", 10**9) < rank.get(x, -1)
                   for x in ("kojic", "arbutin", "hydroquinone") if x in rank)

    c1 = bool(np.isfinite(sc_rho) and sc_rho >= 0.30)
    c2 = bool(np.isfinite(ef) and ef >= 2.0)

    cov = n_sc / max(len(scored), 1)
    print(f"n scored         : {len(scored)} / {len(bench)}   (size-corr. on n={n_sc}, "
          f"coverage {cov:.0%})")
    if cov < 0.80:
        print(f"    !! LOW COVERAGE — rho is over a self-selected subset of {n_sc} "
              f"ligands, not comparable to full-set rho")
    print(f"1'. size-corr rho: {sc_rho:.3f} (p={sc_p:.1e})   {'PASS' if c1 else 'FAIL'}")
    print(f"    raw rho      : {raw_rho:.3f} (p={raw_p:.1e})   [reported, not gating]")
    print(f"2.  EF@20%       : {ef:.2f}                {'PASS' if c2 else 'FAIL'}")
    print(f"3.  anchor order : {rank}  {'PASS' if ordering else 'FAIL'}")

    if is_decoy.any():
        dsc = d.loc[is_decoy, "_s"].iloc[0]
        dr = int((ranked._s > dsc).sum()) + 1
        c4 = dr > 0.25 * (len(ranked) + 1)
        print(f"4.  decoy rank   : {dr}/{len(ranked)+1} (outside top 25%)  {'PASS' if c4 else 'FAIL'}")
        verdict = c1 and c2 and ordering and c4
    else:
        if decoy_expected:
            print("4.  decoy rank   : SKIPPED — decoy absent from scored set")
            print("\nGATE: INCOMPLETE — condition 4 untested, verdict withheld")
            return None
        print("4.  decoy rank   : N/A — attempt predates the condition-4 amendment")
        verdict = c1 and c2 and ordering
        print(f"\nGATE: {'PASS (conditions 1-3)' if verdict else 'FAIL — fix the model, do not generate'}")
        return verdict

    print(f"\nGATE: {'PASS — proceed to generative design' if verdict else 'FAIL — fix the model, do not generate'}")
    return verdict
