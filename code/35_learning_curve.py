#!/usr/bin/env python3
"""How many compounds do we need? Learning curve from the mushroom data."""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
RDLogger.DisableLog("rdApp.*")

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
DESC = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
        "NumRotatableBonds", "NumAromaticRings", "FractionCSP3",
        "HeavyAtomCount", "RingCount", "NHOHCount", "NOCount"]
fn = {n: f for n, f in Descriptors.descList if n in DESC}

df = pd.read_csv(f"{CUR}/multitask_train.csv")
d = df[df.task == "mushroom"].reset_index(drop=True)
ms = [Chem.MolFromSmiles(s) for s in d.smiles]
ok = [i for i, m in enumerate(ms) if m is not None]
d = d.iloc[ok].reset_index(drop=True)
X = np.array([[fn[n](ms[i]) for n in DESC] for i in ok])
y = d.pAct.values
heavy = X[:, DESC.index("HeavyAtomCount")]
print(f"mushroom pool: {len(y)} compounds")

def score(Xs, ys, hs, seed):
    oof = np.zeros(len(ys))
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(Xs):
        m = RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1).fit(Xs[tr], ys[tr])
        oof[te] = m.predict(Xs[te])
    raw = spearmanr(oof, ys).statistic
    sl = linregress(hs, oof)
    corr = spearmanr(oof - (sl.slope * hs + sl.intercept), ys).statistic
    return raw, corr

SIZES = [40, 60, 100, 150, 250, 400, 650, 1000, 1500]
REPS = 8
rng = np.random.default_rng(42)
print(f"\n{'n':>6} {'raw rho':>16} {'size-corrected rho':>22}")
rows = []
for n in SIZES:
    if n > len(y): continue
    raws, corrs = [], []
    for rep in range(REPS):
        idx = rng.choice(len(y), n, replace=False)
        r, c = score(X[idx], y[idx], heavy[idx], 42 + rep)
        raws.append(r); corrs.append(c)
    rows.append({"n": n, "raw_mean": np.mean(raws), "raw_sd": np.std(raws),
                 "corr_mean": np.mean(corrs), "corr_sd": np.std(corrs)})
    print(f"{n:6d}   {np.mean(raws):+.3f} +/- {np.std(raws):.3f}"
          f"      {np.mean(corrs):+.3f} +/- {np.std(corrs):.3f}")

lc = pd.DataFrame(rows)
lc.to_csv(f"{CUR}/learning_curve.csv", index=False)

print("\n== extrapolation: rho vs log(n) ==")
sl = linregress(np.log10(lc.n), lc.corr_mean)
print(f"  size-corrected rho ~ {sl.slope:.3f} * log10(n) + {sl.intercept:.3f}"
      f"   (R^2 = {sl.rvalue**2:.2f})")
for target in (0.55, 0.60, 0.65):
    need = 10 ** ((target - sl.intercept) / sl.slope)
    print(f"  to reach size-corrected rho = {target:.2f}: ~{need:,.0f} compounds")

print(f"\n  your human enzyme set today: 48 compounds, size-corrected rho 0.495")
print("\nCaveat: mushroom and human chemical spaces differ, so this estimates")
print("the shape of the human curve, it does not guarantee it. Treat the")
print("compound counts as an order of magnitude, not a budget line.")
