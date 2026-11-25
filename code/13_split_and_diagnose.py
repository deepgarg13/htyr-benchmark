#!/usr/bin/env python3
"""Split the benchmark into dev/holdout, then diagnose the gate failure on DEV ONLY."""
import os, glob, hashlib, datetime
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
CUR, RES, MOD = f"{ROOT}/data/curated", f"{ROOT}/results", f"{ROOT}/models"
SEED = 20261123

bench = pd.read_csv(f"{CUR}/benchmark_v2.csv")
scored = bench[bench.pAct_median.notna()].copy()

if not os.path.exists(f"{CUR}/dev.csv"):
    rng = np.random.default_rng(SEED)
    idx = rng.permutation(len(scored))
    n_dev = int(round(0.65 * len(scored)))
    dev = scored.iloc[idx[:n_dev]].copy()
    hold = scored.iloc[idx[n_dev:]].copy()
    decoy = bench[bench.pAct_median.isna()]
    dev = pd.concat([dev, decoy], ignore_index=True)   # decoy stays in dev
    dev.to_csv(f"{CUR}/dev.csv", index=False)
    hold.to_csv(f"{CUR}/holdout.csv", index=False)
    sha = hashlib.sha256(open(f"{CUR}/holdout.csv", 'rb').read()).hexdigest()
    with open(f"{ROOT}/BENCHMARK_LOCKED.md", "a") as f:
        f.write(f"""

## AMENDMENT — {datetime.date.today().strftime('%d-%m-%Y')} (dev/holdout split)
After the first gate attempt FAILED on all three conditions, the scored set was
split to protect against overfitting during model iteration.
Random split, seed {SEED}, declared before inspecting membership.
  dev     : {len(dev)} compounds (incl. the unscored decoy) — iterate freely
  holdout : {len(hold)} compounds — SEALED
  holdout SHA256: {sha}

Rule: the holdout is opened ONCE, for the first approach that passes the gate on
dev. Gate criteria are unchanged. Every attempt is logged in ATTEMPTS.md,
pass or fail, before its result is known to be good.
""")
    print(f"split: dev {len(dev)}  holdout {len(hold)} (sealed)")
    if not os.path.exists(f"{ROOT}/ATTEMPTS.md"):
        open(f"{ROOT}/ATTEMPTS.md", "w").write(
            "# Gate attempts\n\n"
            "| # | date | approach | rho (dev) | EF@20% | anchors | decoy | verdict |\n"
            "|---|------|----------|-----------|--------|---------|-------|---------|\n"
            "| 1 | 23-11-2026 | smina vina / ad4 q0 / ad4 q2, neutral ligands, "
            "24A box, deoxy Cu-Cu 4.36A | 0.161 / 0.086 / 0.101 | 1.32 / 1.10 / 1.10 "
            "| arbutin over-ranked | - | FAIL |\n")
        print("created ATTEMPTS.md")

dev = pd.read_csv(f"{CUR}/dev.csv")
dev_ids = set(dev.molecule_chembl_id)
print(f"\ndiagnostics on DEV only ({len(dev)} compounds)\n" + "=" * 58)

def _heavy(s):
    if not isinstance(s, str) or not s.strip():
        return np.nan
    m = Chem.MolFromSmiles(s)
    return m.GetNumHeavyAtoms() if m is not None else np.nan

dev["heavy"] = dev.smiles.map(_heavy)
print(f"  heavy-atom counts: {int(dev.heavy.notna().sum())}/{len(dev)}"
      f"   rows with no SMILES: {int(dev.smiles.isna().sum())}")
for _, r in dev[dev.heavy.isna()].iterrows():
    print(f"    no structure: {r.molecule_chembl_id}  {r.pref_name}")

print("\n== D1: is the score just measuring size? ==")
for run in ("A_vina", "B_ad4_q0", "C_ad4_q2"):
    f = f"{RES}/scores_{run}.csv"
    if not os.path.exists(f): continue
    d = dev.merge(pd.read_csv(f), on="molecule_chembl_id")
    ds = d.dropna(subset=["score", "heavy"])
    dp = d.dropna(subset=["score", "pAct_median"])
    r_size, p_size = spearmanr(ds.score, ds.heavy)
    r_pot, _ = spearmanr(-dp.score, dp.pAct_median)
    print(f"  {run:10s} score~heavy rho = {r_size:+.3f} (p={p_size:.1e}, n={len(ds)})"
          f"   |   score~potency rho = {r_pot:+.3f} (n={len(dp)})")

print("\n== D2: where are the poses actually landing? ==")
cu = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
               for l in open(f"{MOD}/hTYR_receptor.pdbqt")
               if l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu")])

def first_pose(path):
    xyz, on = [], False
    for l in open(path):
        if l.startswith("MODEL"): on = True
        elif l.startswith("ENDMDL") and xyz: break
        elif on and l.startswith(("ATOM", "HETATM")):
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(xyz) if xyz else None

near = 0; total = 0
WATCH = {"CHEMBL450195": "rucinol", "CHEMBL287556": "kojic", "CHEMBL232202": "arbutin",
         "CHEMBL537": "hydroquinone", "CHEMBL43065": "oxyresveratrol"}
for mid in dev.molecule_chembl_id:
    tag = mid.replace(":", "_")
    p = f"{RES}/A_vina_{tag}_out.pdbqt"
    if not os.path.exists(p): continue
    xyz = first_pose(p)
    if xyz is None: continue
    dmin = min(np.linalg.norm(xyz - c, axis=1).min() for c in cu)
    total += 1
    if dmin < 3.0: near += 1
    if mid in WATCH:
        print(f"  {WATCH[mid]:14s} closest atom to Cu: {dmin:5.2f} A")
    if "Thiamidol" in str(dev.loc[dev.molecule_chembl_id == mid, "pref_name"].values):
        print(f"  {'thiamidol':14s} closest atom to Cu: {dmin:5.2f} A")

print(f"\n  poses with any atom within 3.0 A of a copper: {near}/{total} "
      f"({100*near/max(total,1):.0f}%)")
print("  (if this is low, the 24 A box is finding surface pockets, not the active site)")
