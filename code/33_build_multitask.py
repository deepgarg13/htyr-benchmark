#!/usr/bin/env python3
"""Workstream 1, step 1: assemble the multi-fidelity dataset.

Four tasks from files already on disk. The one rule that matters: no compound
in the sealed holdout may appear in ANY task's training rows, or the multi-task
setup leaks and the holdout stops meaning anything.
"""
import os, json
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from chembl_webresource_client.new_client import new_client
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
RAW, CUR = f"{ROOT}/data/raw", f"{ROOT}/data/curated"
CACHE = f"{CUR}/_smiles_cache.json"
FUNGAL = {"CHEMBL3318", "CHEMBL3988582", "CHEMBL6066154"}
POTENCY = {"IC50", "Ki", "EC50", "Kd"}

def ikey(smi):
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) and smi.strip() else None
    if m is None: return None
    try:    return Chem.InchiToInchiKey(Chem.MolToInchi(m)).split("-")[0]
    except Exception: return None

cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
def fetch_smiles(ids):
    want = [i for i in dict.fromkeys(ids) if i.startswith("CHEMBL") and i not in cache]
    if not want: return
    print(f"  fetching {len(want)} structures from ChEMBL", flush=True)
    mol = new_client.molecule
    for i in range(0, len(want), 50):
        for m in mol.filter(molecule_chembl_id__in=want[i:i+50]).only(
                ["molecule_chembl_id", "molecule_structures"]):
            s = (m.get("molecule_structures") or {})
            cache[m["molecule_chembl_id"]] = s.get("canonical_smiles")
        print(f"    {min(i+50, len(want))}/{len(want)}", flush=True)
    json.dump(cache, open(CACHE, "w"))

b = pd.read_csv(f"{CUR}/benchmark_v2.csv")
t1 = b[b.pAct_median.notna() & b.smiles.notna()][
        ["molecule_chembl_id", "smiles", "pAct_median"]].copy()
t1.columns = ["cid", "smiles", "pAct"]; t1["task"] = "human_enzyme"

tri = pd.read_csv(f"{CUR}/human_tyr_triaged.csv")
allact = pd.read_csv(f"{RAW}/tyrosinase_activities.csv")

def from_activities(df, label):
    d = df[df.standard_type.isin(POTENCY)
           & (df.standard_relation == "=")
           & (df.standard_units == "nM")
           & (df.standard_value > 0)].copy()
    d["pAct"] = d.pchembl_value.fillna(-np.log10(d.standard_value * 1e-9))
    g = d.groupby("molecule_chembl_id").pAct.median().reset_index()
    g.columns = ["cid", "pAct"]; g["task"] = label
    return g

t2 = from_activities(tri[tri.bucket == "human_cell"], "human_cell")
t4 = from_activities(allact[allact.target_chembl_id.isin(FUNGAL)], "mushroom")

print("task sizes before structures and de-leaking")
for name, t in (("human_enzyme", t1), ("human_cell", t2), ("mushroom", t4)):
    print(f"  {name:14s} {len(t):5d}")

fetch_smiles(list(t2.cid) + list(t4.cid))
for t in (t2, t4):
    t["smiles"] = t.cid.map(cache)

full = pd.concat([t1, t2[t1.columns], t4[t1.columns]], ignore_index=True)
full = full[full.smiles.notna()].copy()
full["ikey"] = full.smiles.map(ikey)
full = full[full.ikey.notna()]

hold = pd.read_csv(f"{CUR}/holdout.csv")
hold_keys = set(hold.smiles.map(ikey).dropna())
print(f"\nholdout compounds: {len(hold_keys)} unique structures")

leaks = full[full.ikey.isin(hold_keys)].groupby("task").size()
print("holdout structures found in each task (all removed from training):")
for task in ("human_enzyme", "human_cell", "mushroom"):
    print(f"  {task:14s} {int(leaks.get(task, 0)):5d}")

train = full[~full.ikey.isin(hold_keys)].drop_duplicates(["ikey", "task"]).copy()

print(f"\n{'task':16s} {'compounds':>10s} {'pAct range':>18s}")
for task, g in train.groupby("task"):
    print(f"  {task:14s} {len(g):10d} {g.pAct.min():8.2f} - {g.pAct.max():5.2f}")
print(f"  {'TOTAL':14s} {len(train):10d}   ({train.ikey.nunique()} unique structures)")

overlap = (train.groupby("ikey").task.nunique() > 1).sum()
print(f"\ncompounds measured in more than one task: {overlap}")
paired = train[train.task.isin(["human_enzyme", "mushroom"])]
both = paired.groupby("ikey").task.nunique()
print(f"  with BOTH human enzyme and mushroom values: {int((both > 1).sum())}"
      "   <- the training set for Workstream 2")

train.to_csv(f"{CUR}/multitask_train.csv", index=False)
print(f"\nwrote {CUR}/multitask_train.csv")
print("Holdout untouched. Train on this; evaluate only with gate.py.")
