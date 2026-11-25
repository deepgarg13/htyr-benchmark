#!/usr/bin/env python3
"""How much of BindingDB's P14679 data is NEW relative to our locked set?"""
import os, json, re
import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
CUR, RAW = f"{ROOT}/data/curated", f"{ROOT}/data/raw"

data = json.load(open(f"{CUR}/_bindingdb_raw.txt"))
affs = next(iter(data.values())).get("affinities", [])
print(f"BindingDB records: {len(affs)}")
if affs:
    print(f"fields: {sorted(affs[0].keys())}")

df = pd.DataFrame(affs)
for c in ("query", "affinity_type"):
    if c in df:
        print(f"\n{c} values:")
        print(df[c].value_counts().head(10).to_string())

def key(smi):
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    return Chem.InchiToInchiKey(Chem.MolToInchi(m)).split("-")[0] if m else None

df["val"] = pd.to_numeric(df.affinity.astype(str).str.replace(r"[<>=~]", "", regex=True),
                          errors="coerce")
df["censored"] = df.affinity.astype(str).str.contains(r"[<>]")
df["ikey"] = df.smile.map(key)
good = df[df.ikey.notna() & df.val.notna() & ~df.censored]
print(f"\nparsed, exact values, valid structures: {len(good)} records, "
      f"{good.ikey.nunique()} unique compounds")

bench = pd.read_csv(f"{CUR}/benchmark_v2.csv")
bench["ikey"] = bench.smiles.map(key)
have = set(bench.ikey.dropna())
new = good[~good.ikey.isin(have)]
print(f"  already in the locked 74+anchors : {good.ikey.isin(have).sum()} records")
print(f"  NEW compounds                    : {new.ikey.nunique()}")

if "pmid" in good:
    allact = pd.read_csv(f"{RAW}/tyrosinase_activities.csv")
    print(f"\n  distinct source papers in BindingDB : {good.pmid.nunique()}")
    print(f"  distinct papers among NEW compounds : {new.pmid.nunique()}")
    print("\n  top papers contributing new compounds:")
    for pmid, n in new.groupby("pmid").ikey.nunique().nlargest(8).items():
        doi = new[new.pmid == pmid].doi.iloc[0]
        print(f"    PMID {pmid}  {n:3d} compounds   {doi}")

out = new.drop_duplicates("ikey")[["monomerid", "smile", "affinity_type",
                                   "affinity", "val", "pmid", "doi", "ikey"]]
out.to_csv(f"{CUR}/bindingdb_new.csv", index=False)
print(f"\nwrote {CUR}/bindingdb_new.csv ({len(out)} compounds)")
print("\nNOT merged. BindingDB maps by UniProt but does not guarantee the assay")
print("used human enzyme — the source papers still need checking before use.")
