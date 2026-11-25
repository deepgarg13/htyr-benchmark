#!/usr/bin/env python3
"""Ligand-based ranking baseline. Cross-validated on DEV only.
Holdout stays sealed until exactly one model is chosen."""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.model_selection import RepeatedKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
CUR, RES = f"{ROOT}/data/curated", f"{ROOT}/results"
dev = pd.read_csv(f"{CUR}/dev.csv")
dev = dev[dev.smiles.notna() & dev.pAct_median.notna()].reset_index(drop=True)
mols = [Chem.MolFromSmiles(s) for s in dev.smiles]
keep = [i for i, m in enumerate(mols) if m is not None]
dev, mols = dev.iloc[keep].reset_index(drop=True), [mols[i] for i in keep]
y = dev.pAct_median.values
print(f"dev: {len(dev)} compounds, pAct {y.min():.2f}-{y.max():.2f}")

DESC = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
        "NumRotatableBonds", "NumAromaticRings", "FractionCSP3",
        "HeavyAtomCount", "RingCount", "NHOHCount", "NOCount"]
fn = {n: f for n, f in Descriptors.descList if n in DESC}
X_desc = np.array([[fn[n](m) for n in DESC] for m in mols])
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
fps = [gen.GetFingerprint(m) for m in mols]
X_fp = np.array([list(fp) for fp in fps], dtype=float)
heavy = X_desc[:, DESC.index("HeavyAtomCount")]

def evaluate(pred, label):
    r, p = spearmanr(pred, y)
    from scipy.stats import linregress
    sl = linregress(heavy, pred)
    rc, pc = spearmanr(pred - (sl.slope * heavy + sl.intercept), y)
    k = max(1, int(round(0.20 * len(y))))
    top = np.argsort(-pred)[:k]
    ef = (y[top] >= 5.5).mean() / max((y >= 5.5).mean(), 1e-9)
    print(f"  {label:28s} rho={r:+.3f} (p={p:.1e})  size-corr={rc:+.3f} (p={pc:.1e})  EF={ef:.2f}")
    return rc

print("\n== reference points ==")
evaluate(heavy.astype(float), "NULL: heavy atoms only")

MODELS = {
 "Ridge / descriptors":  (make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 3, 20))), X_desc),
 "RF / descriptors":     (RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=42, n_jobs=-1), X_desc),
 "RF / Morgan fp":       (RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=42, n_jobs=-1), X_fp),
 "RF / desc + fp":       (RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=42, n_jobs=-1), np.hstack([X_desc, X_fp])),
}

print("\n== out-of-fold CV on dev (5-fold x 5 repeats) ==")
best = (None, -9)
for name, (model, X) in MODELS.items():
    oof = np.zeros((len(y), 0))
    preds = []
    for rep, (tr, te) in enumerate(RepeatedKFold(n_splits=5, n_repeats=5, random_state=42).split(X)):
        pass
    acc = np.zeros(len(y)); cnt = np.zeros(len(y))
    for tr, te in RepeatedKFold(n_splits=5, n_repeats=5, random_state=42).split(X):
        model.fit(X[tr], y[tr])
        acc[te] += model.predict(X[te]); cnt[te] += 1
    rc = evaluate(acc / cnt, name)
    if rc > best[1]: best = (name, rc)

# Tanimoto k-NN
print()
sim = np.array([DataStructs.BulkTanimotoSimilarity(f, fps) for f in fps])
np.fill_diagonal(sim, -1)
for k in (3, 5):
    pred = np.array([y[np.argsort(-sim[i])[:k]].mean() for i in range(len(y))])
    rc = evaluate(pred, f"Tanimoto {k}-NN")
    if rc > best[1]: best = (f"Tanimoto {k}-NN", rc)

print(f"\nbest on dev by size-corrected rho: {best[0]} ({best[1]:+.3f})")
print("HOLDOUT NOT TOUCHED. Choose one model, then open it once.")
