#!/usr/bin/env python3
"""Workstream 3: interaction fingerprints from the validated receptor.
Hand-rolled, so every feature means something we chose."""
import os, re
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD, RES, CUR = f"{ROOT}/models", f"{ROOT}/results", f"{ROOT}/data/curated"

# Cu-binding histidines, then the hsTYR-specific hot spots named by Roulier 2020
KEY = {180: "HIS", 202: "HIS", 211: "HIS", 363: "HIS", 367: "HIS", 390: "HIS",
       377: "VAL", 380: "SER", 375: "SER", 364: "ASN", 368: "ILE", 352: "THR",
       304: "HIS", 306: "LYS", 308: "ARG", 347: "PHE", 203: "GLU"}

prot, cu = [], []
for l in open(f"{MOD}/hTYR_receptor.pdb"):
    if not l.startswith(("ATOM", "HETATM")): continue
    name, resn, rs = l[12:16].strip(), l[17:20].strip(), l[22:26].strip()
    xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    if resn == "CU": cu.append(xyz); continue
    if name.startswith("H") or l[76:78].strip() == "H": continue
    prot.append((resn, int(rs) if rs.lstrip("-").isdigit() else -1, name[0], xyz))
cu = np.array(cu)[:2]
mid = cu.mean(0)

print("== residue numbering check ==")
found = {r: n for n, r, _, _ in prot}
bad = [(r, KEY[r], found.get(r, "ABSENT")) for r in sorted(KEY) if found.get(r) != KEY[r]]
for r, want, got in bad:
    print(f"  !! residue {r}: expected {want}, found {got}")
if not bad:
    print(f"  all {len(KEY)} key residues match the expected amino acid")
else:
    print("  numbering may be offset - features below will be meaningless")

P = np.array([a[3] for a in prot])
P_res = np.array([a[1] for a in prot])
P_el = np.array([a[2] for a in prot])
KEYS = sorted(KEY)

def pose(path):
    xyz, el, on = [], [], False
    for l in open(path):
        if l.startswith("MODEL"): on = True
        elif l.startswith("ENDMDL") and xyz: break
        elif on and l.startswith(("ATOM", "HETATM")):
            t = l.split()[-1]
            if t == "HD": continue
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            el.append({"A": "C", "OA": "O", "NA": "N", "SA": "S"}.get(t, t[0]))
    return (np.array(xyz), np.array(el)) if xyz else (None, None)

def features(xyz, el):
    f = {}
    for i, tag in enumerate("AB"):
        d = np.linalg.norm(xyz - cu[i], axis=1)
        f[f"dCu{tag}_min"] = d.min()
        for e in "ONS":
            m = el == e
            f[f"dCu{tag}_{e}"] = d[m].min() if m.any() else 12.0
    f["n_near_cu"] = int((np.minimum(np.linalg.norm(xyz - cu[0], axis=1),
                                     np.linalg.norm(xyz - cu[1], axis=1)) < 4.0).sum())
    f["depth"] = np.linalg.norm(xyz.mean(0) - mid)
    D = np.linalg.norm(xyz[:, None, :] - P[None, :, :], axis=2)
    close = D < 4.5
    f["contacts_total"] = int(close.sum())
    polar_l = np.isin(el, ["O", "N"])[:, None]
    polar_p = np.isin(P_el, ["O", "N"])[None, :]
    f["polar_contacts"] = int(((D < 3.5) & polar_l & polar_p).sum())
    for r in KEYS:
        f[f"res_{r}"] = int(close[:, P_res == r].sum())
    return f

df = pd.read_csv(f"{CUR}/multitask_train.csv")
hum = df[df.task == "human_enzyme"].reset_index(drop=True)
rows, keep = [], []
for i, cid in enumerate(hum.cid):
    p = f"{RES}/att4_{cid.replace(':', '_')}_out.pdbqt"
    if not os.path.exists(p): continue
    xyz, el = pose(p)
    if xyz is None: continue
    rows.append(features(xyz, el)); keep.append(i)
hum = hum.iloc[keep].reset_index(drop=True)
IFP = pd.DataFrame(rows)
print(f"\ninteraction fingerprints: {len(IFP)} compounds x {IFP.shape[1]} features")

DESC = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors",
        "NumRotatableBonds", "NumAromaticRings", "FractionCSP3",
        "HeavyAtomCount", "RingCount", "NHOHCount", "NOCount"]
fn = {n: f for n, f in Descriptors.descList if n in DESC}
Xd = np.array([[fn[n](Chem.MolFromSmiles(s)) for n in DESC] for s in hum.smiles])
Xi = IFP.values
Xm = IFP[[c for c in IFP.columns if c.startswith("dCu") or c in
          ("n_near_cu", "depth", "polar_contacts")]].values
y = hum.pAct.values
heavy = Xd[:, DESC.index("HeavyAtomCount")]

def score(X, seed):
    oof = np.zeros(len(y))
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        m = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1).fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    sl = linregress(heavy, oof)
    return spearmanr(oof - (sl.slope * heavy + sl.intercept), y).statistic

print(f"\n== dev CV, n={len(y)}, 20 seeds, matched protocol ==")
for label, X in (("A: descriptors only", Xd),
                 ("B: descriptors + full IFP", np.hstack([Xd, Xi])),
                 ("C: full IFP alone", Xi),
                 ("D: descriptors + metal features", np.hstack([Xd, Xm])),
                 ("E: metal features alone", Xm)):
    s = [score(X, 42 + r) for r in range(20)]
    print(f"  {label:34s} rho = {np.mean(s):+.3f} +/- {np.std(s):.3f}")

m = RandomForestRegressor(n_estimators=600, min_samples_leaf=2,
                          random_state=42, n_jobs=-1).fit(Xi, y)
print("\n== which interactions carry signal? ==")
for n, imp in sorted(zip(IFP.columns, m.feature_importances_),
                     key=lambda t: -t[1])[:12]:
    tag = "  <- hsTYR-specific" if n.startswith("res_") and int(n[4:]) in (
        377, 380, 375, 364, 368, 352, 304, 306, 308) else ""
    print(f"  {n:18s} {imp:.3f}{tag}")
