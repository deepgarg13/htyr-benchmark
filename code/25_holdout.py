#!/usr/bin/env python3
"""THE SINGLE HOLDOUT EVALUATION. One model, one shot, recorded either way."""
import os, datetime
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
CUR = f"{ROOT}/data/curated"
DESC = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
        "NumRotatableBonds", "NumAromaticRings", "FractionCSP3",
        "HeavyAtomCount", "RingCount", "NHOHCount", "NOCount"]
fn = {n: f for n, f in Descriptors.descList if n in DESC}

def featurise(df):
    df = df[df.smiles.notna()].copy()
    ms = [Chem.MolFromSmiles(s) for s in df.smiles]
    ok = [i for i, m in enumerate(ms) if m is not None]
    df = df.iloc[ok].reset_index(drop=True); ms = [ms[i] for i in ok]
    X = np.array([[fn[n](m) for n in DESC] for m in ms])
    return df, X

dev, Xd = featurise(pd.read_csv(f"{CUR}/dev.csv"))
dev_s = dev[dev.pAct_median.notna()]
Xd_s = Xd[dev.pAct_median.notna().values]
hold, Xh = featurise(pd.read_csv(f"{CUR}/holdout.csv"))

model = RandomForestRegressor(n_estimators=500, min_samples_leaf=2,
                              random_state=42, n_jobs=-1)
model.fit(Xd_s, dev_s.pAct_median.values)

print("== what did it learn? (feature importance) ==")
for n, imp in sorted(zip(DESC, model.feature_importances_), key=lambda t: -t[1]):
    print(f"  {n:22s} {imp:.3f}")

pred = model.predict(Xh)
y = hold.pAct_median.values
heavy = Xh[:, DESC.index("HeavyAtomCount")]

raw, p_raw = spearmanr(pred, y)
sl = linregress(heavy, pred)
resid = pred - (sl.slope * heavy + sl.intercept)
corr, p_corr = spearmanr(resid, y)
k = max(1, int(round(0.20 * len(y))))
rate = (y >= 5.5).mean()
ef = (y[np.argsort(-pred)[:k]] >= 5.5).mean() / max(rate, 1e-9)

rng = np.random.default_rng(42)
perm = np.array([spearmanr(resid, rng.permutation(y)).statistic for _ in range(10000)])
p_perm = float((np.abs(perm) >= abs(corr)).mean())

print(f"\n{'='*60}\nHOLDOUT (n={len(y)}, never seen, never tuned on)\n{'='*60}")
print(f"  raw rho              {raw:+.3f} (p={p_raw:.1e})")
print(f"  size-corrected rho   {corr:+.3f} (p={p_corr:.1e}, permutation p={p_perm:.4f})")
print(f"     gate condition 1' >= 0.30   {'PASS' if corr >= 0.30 else 'FAIL'}")
print(f"  EF@20% (top {k} of {len(y)})   {ef:.2f}   active rate {rate:.2f}")
print(f"     gate condition 2  >= 2.0    {'PASS' if ef >= 2.0 else 'FAIL'}")
verdict = corr >= 0.30 and ef >= 2.0
print(f"\n  VERDICT: {'PASS' if verdict else 'FAIL'}")

with open(f"{ROOT}/ATTEMPTS.md", "a") as f:
    f.write(f"| 5 | {datetime.date.today().strftime('%d-%m-%Y')} | RF on 12 RDKit "
            f"descriptors, trained on dev, SINGLE holdout evaluation | "
            f"{corr:+.3f} (holdout) | {ef:.2f} | - | - | "
            f"{'PASS' if verdict else 'FAIL'} |\n")
print(f"\nrecorded in ATTEMPTS.md")
