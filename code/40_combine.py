#!/usr/bin/env python3
"""Structure carries signal but dilutes the model. Combine it properly:
nested feature selection, and prediction-level ensembling."""
import os, importlib.util
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

spec = importlib.util.spec_from_file_location(
    "fx", os.path.expanduser("~/hTYR_Platform/code/39_fix_numbering.py"))
print("reusing feature extraction from script 39...\n")
# rerun 39's extraction by import is messy - rebuild the three matrices instead
exec(open(os.path.expanduser("~/hTYR_Platform/code/39_fix_numbering.py")).read()
     .split("def score(X, seed):")[0])

def oof_pred(X, seed, k=None, Xsel=None):
    """k: if set, pick top-k columns of Xsel by importance on TRAIN only."""
    out = np.zeros(len(y))
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        if k is None:
            Xt, Xe = X[tr], X[te]
        else:
            imp = RandomForestRegressor(n_estimators=200, min_samples_leaf=2,
                        random_state=seed, n_jobs=-1).fit(Xsel[tr], y[tr]).feature_importances_
            cols = np.argsort(-imp)[:k]
            Xt = np.hstack([X[tr], Xsel[tr][:, cols]])
            Xe = np.hstack([X[te], Xsel[te][:, cols]])
        m = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1).fit(Xt, y[tr])
        out[te] = m.predict(Xe)
    return out

def rho(pred):
    sl = linregress(heavy, pred)
    return spearmanr(pred - (sl.slope * heavy + sl.intercept), y).statistic

SEEDS = range(42, 62)
print(f"== dev CV, n={len(y)}, 20 seeds ==")
base = [rho(oof_pred(Xd, s)) for s in SEEDS]
print(f"  A: descriptors only                  rho = {np.mean(base):+.3f} +/- {np.std(base):.3f}")

for k in (1, 2, 3, 5):
    sc = [rho(oof_pred(Xd, s, k=k, Xsel=Xr)) for s in SEEDS]
    print(f"  B{k}: descriptors + top-{k} residue (nested) rho = {np.mean(sc):+.3f} +/- {np.std(sc):.3f}")

print()
for w in (0.2, 0.3, 0.5):
    sc = []
    for s in SEEDS:
        pd_, ps = oof_pred(Xd, s), oof_pred(Xr, s)
        z = lambda v: (v - v.mean()) / (v.std() + 1e-9)
        sc.append(rho((1 - w) * z(pd_) + w * z(ps)))
    print(f"  E{w}: ensemble, {int(100*w)}% structure      rho = {np.mean(sc):+.3f} +/- {np.std(sc):.3f}")

print("\nNested selection means the features are chosen inside each training fold,")
print("so the comparison is honest rather than flattered by seeing all the data.")
