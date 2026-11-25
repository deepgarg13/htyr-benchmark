#!/usr/bin/env python3
"""Put the human set on the same protocol as the curve, and ask whether it is
an easier problem or was just measured more kindly."""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
RDLogger.DisableLog("rdApp.*")

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
DESC = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
        "NumRotatableBonds", "NumAromaticRings", "FractionCSP3",
        "HeavyAtomCount", "RingCount", "NHOHCount", "NOCount"]
fn = {n: f for n, f in Descriptors.descList if n in DESC}
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

def prep(d):
    ms = [Chem.MolFromSmiles(s) for s in d.smiles]
    ok = [i for i, m in enumerate(ms) if m is not None]
    d = d.iloc[ok].reset_index(drop=True); ms = [ms[i] for i in ok]
    X = np.array([[fn[n](m) for n in DESC] for m in ms])
    return d, X, [gen.GetFingerprint(m) for m in ms]

df = pd.read_csv(f"{CUR}/multitask_train.csv")
hum, Xh, fh = prep(df[df.task == "human_enzyme"])
mush, Xm, fm = prep(df[df.task == "mushroom"])

def score_once(X, y, h, seed):
    oof = np.zeros(len(y))
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        m = RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1).fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    sl = linregress(h, oof)
    return spearmanr(oof - (sl.slope * h + sl.intercept), y).statistic

hv = Xh[:, DESC.index("HeavyAtomCount")]
hs = [score_once(Xh, hum.pAct.values, hv, 42 + r) for r in range(20)]
print(f"human n={len(hum)}, SAME protocol as the curve")
print(f"  size-corrected rho = {np.mean(hs):+.3f} +/- {np.std(hs):.3f}"
      f"   (range {min(hs):+.3f} to {max(hs):+.3f})")
print(f"  previously reported (averaged predictions): +0.495")

lc = pd.read_csv(f"{CUR}/learning_curve.csv")
sl = linregress(np.log10(lc.n), lc.corr_mean)
equiv = 10 ** ((np.mean(hs) - sl.intercept) / sl.slope)
print(f"\n  equivalent point on the mushroom curve: n ~ {equiv:,.0f} compounds")

def diversity(fps, idx=None):
    f = [fps[i] for i in idx] if idx is not None else fps
    s = [np.mean(DataStructs.BulkTanimotoSimilarity(f[i], f)) for i in range(len(f))]
    return 1 - np.mean(s)

rng = np.random.default_rng(42)
sub = rng.choice(len(mush), len(hum), replace=False)
print(f"\n== is the human set simply an easier problem? ==")
print(f"{'':22s} {'pAct sd':>9} {'range':>9} {'diversity':>11}")
print(f"  human (n={len(hum):3d})       {hum.pAct.std():9.2f}"
      f" {hum.pAct.max()-hum.pAct.min():9.2f} {diversity(fh):11.3f}")
print(f"  mushroom slice (n={len(sub)}) {mush.pAct.values[sub].std():9.2f}"
      f" {np.ptp(mush.pAct.values[sub]):9.2f} {diversity(fm, sub):11.3f}")
print(f"  mushroom all ({len(mush)})   {mush.pAct.std():9.2f}"
      f" {np.ptp(mush.pAct.values):9.2f}")

ms_sub = [score_once(Xm[sub], mush.pAct.values[sub],
                     Xm[sub][:, DESC.index("HeavyAtomCount")], 42 + r)
          for r in range(20)]
print(f"\n  mushroom slice, same n, same protocol:"
      f" {np.mean(ms_sub):+.3f} +/- {np.std(ms_sub):.3f}")
print("\nIf the human set has a wider pAct spread and higher diversity, it is")
print("genuinely easier and the curve's compound counts are an overestimate.")
