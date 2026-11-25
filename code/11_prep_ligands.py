#!/usr/bin/env python3
"""SMILES -> 3D -> pdbqt for the locked benchmark."""
import os, subprocess
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
LIG = os.path.expanduser("~/hTYR_Platform/ligands")
os.makedirs(LIG, exist_ok=True)

df = pd.read_csv(f"{CUR}/benchmark_v2.csv")
df = df[df.smiles.notna()]
print(f"{len(df)} ligands with structures")

ok, bad = 0, []
for _, r in df.iterrows():
    tag = r.molecule_chembl_id.replace(":", "_")
    out = f"{LIG}/{tag}.pdbqt"
    if os.path.exists(out):
        ok += 1; continue
    m = Chem.MolFromSmiles(r.smiles)
    if m is None:
        bad.append((tag, "unparseable SMILES")); continue
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) != 0:
        bad.append((tag, "embedding failed")); continue
    AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    sdf = f"{LIG}/{tag}.sdf"
    Chem.SDWriter(sdf).write(m)
    subprocess.run(["obabel", sdf, "-O", out], capture_output=True)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        ok += 1
    else:
        bad.append((tag, "pdbqt conversion failed"))

print(f"  prepared: {ok}")
for t, why in bad:
    print(f"  FAILED {t}: {why}")
