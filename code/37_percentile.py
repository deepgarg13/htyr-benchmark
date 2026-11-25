#!/usr/bin/env python3
"""Where does the human set sit among equally-sized mushroom subsamples?
And is assay heterogeneity the reason?"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
RAW, CUR = f"{ROOT}/data/raw", f"{ROOT}/data/curated"
DESC = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
        "NumRotatableBonds", "NumAromaticRings", "FractionCSP3",
        "HeavyAtomCount", "RingCount", "NHOHCount", "NOCount"]
fn = {n: f for n, f in Descriptors.descList if n in DESC}

def prep(d):
    ms = [Chem.MolFromSmiles(s) for s in d.smiles]
    ok = [i for i, m in enumerate(ms) if m is not None]
    d = d.iloc[ok].reset_index(drop=True)
    return d, np.array([[fn[n](ms[i]) for n in DESC] for i in ok])

df = pd.read_csv(f"{CUR}/multitask_train.csv")
hum, Xh = prep(df[df.task == "human_enzyme"])
mush, Xm = prep(df[df.task == "mushroom"])
N = len(hum)

def score(X, y, seed):
    oof = np.zeros(len(y))
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        m = RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1).fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    h = X[:, DESC.index("HeavyAtomCount")]
    sl = linregress(h, oof)
    return (spearmanr(oof, y).statistic,
            spearmanr(oof - (sl.slope * h + sl.intercept), y).statistic)

h_raw, h_corr = np.mean([score(Xh, hum.pAct.values, 42 + r) for r in range(10)], axis=0)
print(f"human (n={N}): raw {h_raw:+.3f}   size-corrected {h_corr:+.3f}\n")

DRAWS = 150
rng = np.random.default_rng(7)
raws, corrs = [], []
for i in range(DRAWS):
    idx = rng.choice(len(mush), N, replace=False)
    r, c = score(Xm[idx], mush.pAct.values[idx], 1000 + i)
    raws.append(r); corrs.append(c)
    if (i + 1) % 30 == 0: print(f"  {i+1}/{DRAWS} draws", flush=True)

raws, corrs = np.array(raws), np.array(corrs)
print(f"\nmushroom subsamples at n={N}, {DRAWS} draws")
print(f"  raw            mean {raws.mean():+.3f} +/- {raws.std():.3f}"
      f"   [{np.percentile(raws,5):+.3f}, {np.percentile(raws,95):+.3f}]")
print(f"  size-corrected mean {corrs.mean():+.3f} +/- {corrs.std():.3f}"
      f"   [{np.percentile(corrs,5):+.3f}, {np.percentile(corrs,95):+.3f}]")
print(f"\n  human sits at the {100*(corrs < h_corr).mean():.0f}th percentile"
      f" of mushroom subsamples")
print(f"  beats {int((corrs < h_corr).sum())} of {DRAWS} draws")

print("\n== assay heterogeneity: how many source papers per set? ==")
tri = pd.read_csv(f"{CUR}/human_tyr_triaged.csv")
allact = pd.read_csv(f"{RAW}/tyrosinase_activities.csv")
FUNGAL = {"CHEMBL3318", "CHEMBL3988582", "CHEMBL6066154"}
hdocs = tri[(tri.bucket == "human_enzyme")].groupby("molecule_chembl_id").document_chembl_id.first()
mdocs = (allact[allact.target_chembl_id.isin(FUNGAL)]
         .groupby("molecule_chembl_id").document_chembl_id.first())
nh = hdocs.reindex(hum.cid).dropna().nunique()
print(f"  human set:  {nh} distinct papers for {N} compounds"
      f"  ({N/max(nh,1):.1f} compounds per paper)")
samples = [mdocs.reindex(mush.cid.values[rng.choice(len(mush), N, replace=False)])
           .dropna().nunique() for _ in range(30)]
print(f"  mushroom:   {np.mean(samples):.0f} distinct papers for {N} compounds"
      f"  ({N/max(np.mean(samples),1):.1f} compounds per paper)")
print("\nFewer papers per compound means one protocol, less inter-lab noise,")
print("and congeneric series the model can interpolate within.")
