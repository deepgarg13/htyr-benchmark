"""The frozen gate. Import and call once docking scores exist."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ANCHOR = dict(rucinol="CHEMBL450195", kojic="CHEMBL287556",
              arbutin="CHEMBL232202", hydroquinone="CHEMBL537")

def run_gate(bench, scores, score_col="score", higher_is_better=True):
    """bench: benchmark CSV as DataFrame. scores: molecule_chembl_id -> score."""
    d = bench.merge(scores, on="molecule_chembl_id", how="inner")
    s = d[score_col] if higher_is_better else -d[score_col]
    rho, p = spearmanr(s, d.pAct_median)

    d = d.assign(_s=s).sort_values("_s", ascending=False)
    k = max(1, int(round(0.20 * len(d))))
    act = d.pAct_median >= 5.5
    ef = (d.head(k).pAct_median >= 5.5).sum() / k / (act.mean() or np.nan)

    rank = {n: int((d._s > d.loc[d.molecule_chembl_id == c, "_s"].iloc[0]).sum()) + 1
            for n, c in ANCHOR.items() if (d.molecule_chembl_id == c).any()}
    ordering = all(rank.get("rucinol", 1e9) < rank.get(x, -1)
                   for x in ("kojic", "arbutin", "hydroquinone") if x in rank)

    print(f"n matched        : {len(d)} / {len(bench)}")
    print(f"1. Spearman rho  : {rho:.3f} (p={p:.1e})   {'PASS' if rho >= 0.30 else 'FAIL'}")
    print(f"2. EF@20%        : {ef:.2f}                {'PASS' if ef >= 2.0 else 'FAIL'}")
    print(f"3. anchor order  : {rank}  {'PASS' if ordering else 'FAIL'}")
    decoy = d[d.role == "decoy_mushroom_potent"] if "role" in d else d.iloc[0:0]
    decoy_ok = True
    if not decoy.empty:
        dr = int((d._s > decoy._s.iloc[0]).sum()) + 1
        decoy_ok = dr > 0.25 * len(d)
        print(f"4. decoy rank    : {dr}/{len(d)} (must be outside top 25%)  "
              f"{'PASS' if decoy_ok else 'FAIL'}")
    verdict = (rho >= 0.30) and (ef >= 2.0) and ordering and decoy_ok
    print(f"\nGATE: {'PASS — proceed to generative design' if verdict else 'FAIL — fix the model, do not generate'}")
    return verdict
