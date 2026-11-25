#!/usr/bin/env python3
"""Is the gate passable by a null size model? And does a chemically specific
Cu-coordination filter help? Reuses attempt-2 poses; no re-docking."""
import os, re, sys, glob
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD, RES, CUR = f"{ROOT}/models", f"{ROOT}/results", f"{ROOT}/data/curated"

dev = pd.read_csv(f"{CUR}/dev.csv")
dev["heavy"] = dev.smiles.map(lambda s: Chem.MolFromSmiles(s).GetNumHeavyAtoms()
                              if isinstance(s, str) and Chem.MolFromSmiles(s) else np.nan)
v = dev.dropna(subset=["heavy", "pAct_median"])

print("== THE CONFOUND: does potency track size in this benchmark? ==")
r, p = spearmanr(v.heavy, v.pAct_median)
print(f"  pAct ~ heavy atoms:  rho = {r:+.3f}  (p={p:.1e}, n={len(v)})")

print("\n== NULL MODEL: rank by heavy-atom count alone ==")
r_null, p_null = spearmanr(v.heavy, v.pAct_median)
k = max(1, int(round(0.20 * len(v))))
top = v.nlargest(k, "heavy")
ef_null = (top.pAct_median >= 5.5).mean() / max((v.pAct_median >= 5.5).mean(), 1e-9)
print(f"  size-only rho   = {r_null:+.3f}")
print(f"  size-only EF@20% = {ef_null:.2f}")
if r_null >= 0.30 and ef_null >= 2.0:
    print("  !! A PURE SIZE COUNTER PASSES CONDITIONS 1 AND 2.")
    print("     The gate needs a size-corrected condition before any model is judged by it.")
elif r_null >= 0.30:
    print("  !! a pure size counter passes condition 1 — condition 1 alone is not evidence")
else:
    print("  OK: size alone does not pass. The gate still discriminates.")

print("\n== chemically specific Cu coordination ==")
cu = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
               for l in open(f"{MOD}/hTYR_receptor.pdbqt")
               if l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu")])

def poses(path):
    out, xyz, typ, sc = [], [], [], None
    for l in open(path):
        if l.startswith("MODEL"): xyz, typ, sc = [], [], None
        elif l.startswith("REMARK"):
            m = re.search(r"(?:minimizedAffinity|VINA RESULT:)\s+(-?\d+\.\d+)", l)
            if m: sc = float(m.group(1))
        elif l.startswith(("ATOM", "HETATM")):
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            typ.append(l.split()[-1])
        elif l.startswith("ENDMDL") and xyz and sc is not None:
            out.append((sc, np.array(xyz), typ))
    return out

rows = []
for mid in dev.molecule_chembl_id:
    op = f"{RES}/att2_{mid.replace(':','_')}_out.pdbqt"
    if not os.path.exists(op): continue
    best_o = np.nan
    for sc, xyz, typ in poses(op):
        oxy = [i for i, t in enumerate(typ) if t in ("OA", "O")]
        if not oxy: continue
        dmin = min(np.linalg.norm(xyz[i] - c) for i in oxy for c in cu)
        if dmin < 2.8:
            best_o = sc; break          # poses are score-ordered
    rows.append({"molecule_chembl_id": mid, "score_coord": best_o})

d = dev.merge(pd.DataFrame(rows), on="molecule_chembl_id")
d["heavy"] = d.smiles.map(lambda s: Chem.MolFromSmiles(s).GetNumHeavyAtoms()
                          if isinstance(s, str) and Chem.MolFromSmiles(s) else np.nan)
x = d.dropna(subset=["score_coord", "pAct_median", "heavy"])
print(f"  ligands with an O within 2.8 A of a copper: {len(x)}/{len(d)}")
if len(x) >= 12:
    r_c, p_c = spearmanr(-x.score_coord, x.pAct_median)
    r_s, _ = spearmanr(x.score_coord, x.heavy)
    print(f"  coordinating-pose rho = {r_c:+.3f} (p={p_c:.1e}, n={len(x)})"
          f"   size-rho = {r_s:+.3f}")
else:
    print("  too few coordinating poses to judge")
d.to_csv(f"{RES}/attempt2b_dev.csv", index=False)

print("\n== anchors: does the coordinating filter put them right? ==")
W = {"CHEMBL450195": "rucinol", "CHEMBL287556": "kojic", "CHEMBL232202": "arbutin",
     "CHEMBL537": "hydroquinone", "CHEMBL43065": "oxyresveratrol"}
xs = x.sort_values("score_coord")
for cid, nm in W.items():
    row = xs[xs.molecule_chembl_id == cid]
    if row.empty:
        print(f"  {nm:14s} no coordinating pose")
    else:
        rk = int((xs.score_coord < row.score_coord.iloc[0]).sum()) + 1
        print(f"  {nm:14s} rank {rk}/{len(xs)}   pAct {row.pAct_median.iloc[0]:.2f}")
