#!/usr/bin/env python3
"""How many compounds, across how many series, to score novel scaffolds?
Scaffold-split learning curve on the mushroom pool."""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from sklearn.ensemble import RandomForestRegressor
RDLogger.DisableLog("rdApp.*")

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
DESC = ["MolWt","MolLogP","TPSA","NumHDonors","NumHAcceptors","NumRotatableBonds",
        "NumAromaticRings","FractionCSP3","HeavyAtomCount","RingCount","NHOHCount","NOCount"]
fn = {n: f for n, f in Descriptors.descList if n in DESC}
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

df = pd.read_csv(f"{CUR}/multitask_train.csv")
d = df[df.task == "mushroom"].reset_index(drop=True)
ms = [Chem.MolFromSmiles(s) for s in d.smiles]
ok = [i for i, m in enumerate(ms) if m is not None]
d = d.iloc[ok].reset_index(drop=True); ms = [ms[i] for i in ok]
X = np.array([[fn[n](m) for n in DESC] for m in ms])
y = d.pAct.values
heavy = X[:, DESC.index("HeavyAtomCount")]
fps = [gen.GetFingerprint(m) for m in ms]
print(f"pool: {len(y)} compounds")

print("clustering into series (this takes a minute)...", flush=True)
dists = []
for i in range(1, len(fps)):
    dists.extend(1 - np.array(DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])))
cl = Butina.ClusterData(dists, len(fps), 0.5, isDistData=True)
grp = np.zeros(len(fps), dtype=int)
for gi, c in enumerate(cl):
    for i in c: grp[i] = gi
print(f"  {len(cl)} series, largest {max(len(c) for c in cl)}")

def run(idx, seed):
    g = grp[idx]; uniq = np.unique(g)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(uniq))
    fold = {uniq[o]: i % 5 for i, o in enumerate(order)}
    out = np.full(len(idx), np.nan)
    for f in range(5):
        te = np.array([fold[v] == f for v in g])
        if te.sum() < 3 or (~te).sum() < 10: continue
        m = RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1).fit(X[idx][~te], y[idx][~te])
        out[te] = m.predict(X[idx][te])
    k = ~np.isnan(out)
    if k.sum() < 10: return np.nan, len(uniq)
    sl = linregress(heavy[idx][k], out[k])
    return spearmanr(out[k] - (sl.slope*heavy[idx][k] + sl.intercept), y[idx][k]).statistic, len(uniq)

rng = np.random.default_rng(42)
print(f"\n{'n':>6} {'series':>8} {'scaffold-split rho':>22}")
rows = []
for n in (50, 100, 200, 400, 700, 1100, 1600):
    if n > len(y): continue
    sc, gs = [], []
    for rep in range(6):
        idx = rng.choice(len(y), n, replace=False)
        r, g = run(idx, 100 + rep)
        if not np.isnan(r): sc.append(r); gs.append(g)
    rows.append({"n": n, "series": np.mean(gs), "rho": np.mean(sc), "sd": np.std(sc)})
    print(f"{n:6d} {np.mean(gs):8.0f}   {np.mean(sc):+.3f} +/- {np.std(sc):.3f}")

lc = pd.DataFrame(rows); lc.to_csv(f"{CUR}/scaffold_curve.csv", index=False)
sl = linregress(np.log10(lc.n), lc.rho)
print(f"\n  rho ~ {sl.slope:.3f}*log10(n) + {sl.intercept:.3f}  (R^2={sl.rvalue**2:.2f})")
for t in (0.30, 0.40, 0.50):
    n = 10 ** ((t - sl.intercept) / sl.slope)
    per = lc.series.sum() / lc.n.sum()
    print(f"  out-of-series rho {t:.2f}: ~{n:,.0f} compounds, ~{n*per:,.0f} series")
print(f"\n  your human set today: 48 compounds, 25 series, scaffold-split rho 0.191")
