#!/usr/bin/env python3
"""Attempt 2: tight box + Cu-contact pose selection + size correction. DEV ONLY."""
import os, re, json, subprocess, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD, LIG, RES, CUR = f"{ROOT}/models", f"{ROOT}/ligands", f"{ROOT}/results", f"{ROOT}/data/curated"
sys.path.insert(0, f"{ROOT}/code")

box = json.load(open(f"{MOD}/box.json"))
SIZE = int(sys.argv[1]) if len(sys.argv) > 1 else 18
dev = pd.read_csv(f"{CUR}/dev.csv")
cu = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
               for l in open(f"{MOD}/hTYR_receptor.pdbqt")
               if l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu")])

def poses(path):
    out, cur, sc = [], [], None
    for l in open(path):
        if l.startswith("MODEL"): cur, sc = [], None
        elif l.startswith("REMARK"):
            m = re.search(r"(?:minimizedAffinity|VINA RESULT:)\s+(-?\d+\.\d+)", l)
            if m: sc = float(m.group(1))
        elif l.startswith(("ATOM", "HETATM")):
            cur.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        elif l.startswith("ENDMDL"):
            if cur and sc is not None: out.append((sc, np.array(cur)))
    return out

rows = []
print(f"attempt 2: box {SIZE} A, 20 modes, dev only ({len(dev)} compounds)")
for i, mid in enumerate(dev.molecule_chembl_id, 1):
    tag = mid.replace(":", "_")
    lp, op = f"{LIG}/{tag}.pdbqt", f"{RES}/att2_{tag}_out.pdbqt"
    if not os.path.exists(lp): continue
    if not os.path.exists(op):
        subprocess.run(["smina", "-r", f"{MOD}/hTYR_receptor.pdbqt", "-l", lp,
            "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
            "--center_z", str(box["center_z"]), "--size_x", str(SIZE),
            "--size_y", str(SIZE), "--size_z", str(SIZE),
            "--exhaustiveness", "8", "--num_modes", "20", "--seed", "42",
            "--cpu", "4", "-o", op], capture_output=True)
    ps = poses(op)
    if not ps: continue
    best = ps[0][0]
    contact = [(s, x) for s, x in ps
               if min(np.linalg.norm(x - c, axis=1).min() for c in cu) < 3.0]
    rows.append({"molecule_chembl_id": mid, "score_best": best,
                 "score_cu": contact[0][0] if contact else np.nan,
                 "n_cu_poses": len(contact)})
    if i % 10 == 0: print(f"  {i}/{len(dev)}", flush=True)

d = dev.merge(pd.DataFrame(rows), on="molecule_chembl_id")
d["heavy"] = d.smiles.map(lambda s: Chem.MolFromSmiles(s).GetNumHeavyAtoms()
                          if isinstance(s, str) and Chem.MolFromSmiles(s) else np.nan)
print(f"\nligands with >=1 Cu-contacting pose: {(d.n_cu_poses > 0).sum()}/{len(d)}"
      f"  ({100*(d.n_cu_poses > 0).mean():.0f}%)")

d.to_csv(f"{RES}/attempt2_dev.csv", index=False)
v = d.dropna(subset=["heavy", "pAct_median"])

print("\n== four ranking variants (dev) ==")
for label, col in (("2a best pose, raw", "score_best"), ("2b Cu pose, raw", "score_cu")):
    x = v.dropna(subset=[col])
    if len(x) < 10: print(f"  {label:26s} too few ({len(x)})"); continue
    r_pot, p = spearmanr(-x[col], x.pAct_median)
    r_sz, _ = spearmanr(x[col], x.heavy)
    print(f"  {label:26s} rho={r_pot:+.3f} (p={p:.1e}, n={len(x)})  size-rho={r_sz:+.3f}")
    # size-corrected: residual of score on heavy atoms
    sl = linregress(x.heavy, x[col])
    resid = x[col] - (sl.slope * x.heavy + sl.intercept)
    r_res, p_res = spearmanr(-resid, x.pAct_median)
    r_le, _ = spearmanr(-(x[col] / x.heavy), x.pAct_median)
    print(f"  {'  -> size-corrected':26s} rho={r_res:+.3f} (p={p_res:.1e})"
          f"   ligand-efficiency rho={r_le:+.3f}")
