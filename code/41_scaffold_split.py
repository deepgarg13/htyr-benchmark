#!/usr/bin/env python3
"""The test that matters: can the model score a scaffold it has never seen?"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from sklearn.ensemble import RandomForestRegressor
RDLogger.DisableLog("rdApp.*")

exec(open(os.path.expanduser("~/hTYR_Platform/code/39_fix_numbering.py")).read()
     .split("def score(X, seed):")[0])

mols = [Chem.MolFromSmiles(s) for s in hum.smiles]
scaf = [MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols]
print(f"\n{len(set(scaf))} distinct Murcko scaffolds among {len(mols)} compounds")

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
fps = [gen.GetFingerprint(m) for m in mols]
dists = []
for i in range(1, len(fps)):
    dists.extend(1 - np.array(DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])))
for cutoff in (0.4, 0.5, 0.6):
    cl = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
    print(f"  Tanimoto cutoff {cutoff}: {len(cl)} chemical series"
          f"  (largest {max(len(c) for c in cl)})")

cl = Butina.ClusterData(dists, len(fps), 0.5, isDistData=True)
grp = np.zeros(len(fps), dtype=int)
for gi, c in enumerate(cl):
    for i in c: grp[i] = gi
n_groups = len(cl)
print(f"\nusing {n_groups} series at cutoff 0.5")

def grouped_oof(X, seed, folds=5):
    rng = np.random.default_rng(seed)
    order = rng.permutation(n_groups)
    fold_of = {g: i % folds for i, g in enumerate(order)}
    out = np.full(len(y), np.nan)
    for f in range(folds):
        te = np.array([fold_of[g] == f for g in grp])
        if te.sum() < 3 or (~te).sum() < 10: continue
        m = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1).fit(X[~te], y[~te])
        out[te] = m.predict(X[te])
    return out

def rho(pred):
    ok = ~np.isnan(pred)
    sl = linregress(heavy[ok], pred[ok])
    return spearmanr(pred[ok] - (sl.slope * heavy[ok] + sl.intercept), y[ok]).statistic

Xall = np.hstack([Xd, Xr])
SEEDS = range(42, 72)
print(f"\n== SCAFFOLD-SPLIT CV (series held out whole), n={len(y)}, 30 seeds ==")
for label, X in (("A: descriptors only", Xd),
                 ("C: residue contacts alone", Xr),
                 ("D: descriptors + residue contacts", Xall)):
    s = [rho(grouped_oof(X, sd)) for sd in SEEDS]
    print(f"  {label:36s} rho = {np.mean(s):+.3f} +/- {np.std(s):.3f}")

print("\n== the same models under RANDOM split, for comparison ==")
from sklearn.model_selection import KFold
def rand_oof(X, seed):
    out = np.zeros(len(y))
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        m = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1).fit(X[tr], y[tr])
        out[te] = m.predict(X[te])
    return out
for label, X in (("A: descriptors only", Xd),
                 ("C: residue contacts alone", Xr),
                 ("D: descriptors + residue contacts", Xall)):
    s = [rho(rand_oof(X, sd)) for sd in SEEDS]
    print(f"  {label:36s} rho = {np.mean(s):+.3f} +/- {np.std(s):.3f}")

print("\nThe gap between the two blocks is how much of today's performance was")
print("interpolation within series. If structure holds up better under the")
print("scaffold split, that is the case for Workstream 3.")
