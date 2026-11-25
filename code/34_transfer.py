#!/usr/bin/env python3
"""Does mushroom data transfer? Stacking test. DEV ONLY - holdout stays shut."""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RepeatedKFold, KFold
RDLogger.DisableLog("rdApp.*")

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
DESC = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
        "NumRotatableBonds", "NumAromaticRings", "FractionCSP3",
        "HeavyAtomCount", "RingCount", "NHOHCount", "NOCount"]
fn = {n: f for n, f in Descriptors.descList if n in DESC}

df = pd.read_csv(f"{CUR}/multitask_train.csv")
def feats(d):
    ms = [Chem.MolFromSmiles(s) for s in d.smiles]
    ok = [i for i, m in enumerate(ms) if m is not None]
    d = d.iloc[ok].reset_index(drop=True)
    X = np.array([[fn[n](ms[i]) for n in DESC] for i in ok])
    return d, X

mush, Xm = feats(df[df.task == "mushroom"])
hum,  Xh = feats(df[df.task == "human_enzyme"])
print(f"mushroom {len(mush)}   human_enzyme {len(hum)}")

print("\n== does the mushroom model itself work? ==")
oof = np.zeros(len(mush))
for tr, te in KFold(5, shuffle=True, random_state=42).split(Xm):
    m = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                              random_state=42, n_jobs=-1).fit(Xm[tr], mush.pAct.values[tr])
    oof[te] = m.predict(Xm[te])
r, p = spearmanr(oof, mush.pAct)
print(f"  mushroom OOF rho = {r:+.3f} (p={p:.1e}, n={len(mush)})")

mush_model = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                                   random_state=42, n_jobs=-1).fit(Xm, mush.pAct.values)
hum_mush_pred = mush_model.predict(Xh)
r_direct, p_direct = spearmanr(hum_mush_pred, hum.pAct)
print(f"\n  predicted-mushroom vs human pAct, direct: rho = {r_direct:+.3f} (p={p_direct:.1e})")
print("  (a strongly NEGATIVE or near-zero value is expected and is fine -")
print("   the question is whether it adds signal, not whether it agrees)")

paired = df[df.task.isin(["human_enzyme", "mushroom"])]
both = paired.groupby("ikey").filter(lambda g: g.task.nunique() > 1)
if len(both):
    w = both.pivot_table(index="ikey", columns="task", values="pAct")
    rr, pp = spearmanr(w.human_enzyme, w.mushroom)
    print(f"\n  MEASURED human vs mushroom, {len(w)} paired compounds: rho = {rr:+.3f} (p={pp:.2f})")

y = hum.pAct.values
heavy = Xh[:, DESC.index("HeavyAtomCount")]
def evaluate(pred, label):
    r, p = spearmanr(pred, y)
    sl = linregress(heavy, pred)
    rc, pc = spearmanr(pred - (sl.slope * heavy + sl.intercept), y)
    k = max(1, int(round(0.20 * len(y))))
    ef = (y[np.argsort(-pred)[:k]] >= 5.5).mean() / max((y >= 5.5).mean(), 1e-9)
    print(f"  {label:34s} rho={r:+.3f}  size-corr={rc:+.3f} (p={pc:.1e})  EF={ef:.2f}")

def cv(X):
    acc, cnt = np.zeros(len(y)), np.zeros(len(y))
    for tr, te in RepeatedKFold(n_splits=5, n_repeats=5, random_state=42).split(X):
        m = RandomForestRegressor(n_estimators=500, min_samples_leaf=2,
                                  random_state=42, n_jobs=-1).fit(X[tr], y[tr])
        acc[te] += m.predict(X[te]); cnt[te] += 1
    return acc / cnt

print(f"\n== human_enzyme, dev CV (n={len(y)}) ==")
evaluate(cv(Xh), "A: 12 descriptors (baseline)")
evaluate(cv(np.hstack([Xh, hum_mush_pred[:, None]])), "B: + predicted mushroom pAct")
evaluate(hum_mush_pred, "C: predicted mushroom alone")

print("\nB beating A means transfer is real. C is the control:")
print("if C alone already explains it, the mushroom model is doing the work.")
