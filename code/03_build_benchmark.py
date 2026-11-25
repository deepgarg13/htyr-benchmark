#!/usr/bin/env python3
"""Build the locked human-TYR enzyme benchmark: filter, dedupe, attach SMILES."""
import os
import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
df = pd.read_csv(f"{CUR}/human_tyr_triaged.csv")

e = df[(df.bucket == "human_enzyme")
       & (df.standard_relation == "=")
       & (df.standard_units == "nM")
       & (df.standard_value > 0)].copy()
print(f"human_enzyme, exact, nM, positive : {len(e)} measurements, "
      f"{e.molecule_chembl_id.nunique()} compounds")

e["pAct"] = e.pchembl_value.fillna(-np.log10(e.standard_value * 1e-9))

g = e.groupby("molecule_chembl_id").agg(
    pAct_median=("pAct", "median"),
    pAct_spread=("pAct", lambda s: float(s.max() - s.min())),
    n_measurements=("pAct", "size"),
    types=("standard_type", lambda s: "|".join(sorted(set(s)))),
    docs=("document_chembl_id", lambda s: "|".join(sorted(set(s.dropna())))),
).reset_index()

print(f"after dedupe                      : {len(g)} compounds")
noisy = g[(g.n_measurements > 1) & (g.pAct_spread > 1.0)]
print(f"  of which disagree >1 log unit   : {len(noisy)}  (flagged, not dropped)")
g["flag_noisy"] = g.molecule_chembl_id.isin(noisy.molecule_chembl_id)

print("\nfetching structures...")
mol = new_client.molecule
smiles = {}
ids = g.molecule_chembl_id.tolist()
for i in range(0, len(ids), 50):
    chunk = ids[i:i + 50]
    for m in mol.filter(molecule_chembl_id__in=chunk).only(
            ["molecule_chembl_id", "pref_name", "molecule_structures"]):
        s = (m.get("molecule_structures") or {})
        smiles[m["molecule_chembl_id"]] = (m.get("pref_name"),
                                           s.get("canonical_smiles"))
    print(f"  {min(i+50, len(ids))}/{len(ids)}", flush=True)

g["pref_name"] = g.molecule_chembl_id.map(lambda c: (smiles.get(c) or (None, None))[0])
g["smiles"]    = g.molecule_chembl_id.map(lambda c: (smiles.get(c) or (None, None))[1])
g = g.sort_values("pAct_median", ascending=False)
g.to_csv(f"{CUR}/benchmark_human_enzyme.csv", index=False)
print(f"\nwrote {CUR}/benchmark_human_enzyme.csv")

print("\n== potency range ==")
print(f"  strongest : {g.pAct_median.max():.2f}   weakest: {g.pAct_median.min():.2f}")
print(f"  span      : {g.pAct_median.max() - g.pAct_median.min():.2f} log units")
print(f"  missing SMILES: {g.smiles.isna().sum()}")

ANCHORS = {"CHEMBL450195": "rucinol / 4-n-butylresorcinol",
           "CHEMBL602": "cysteamine", "CHEMBL537": "hydroquinone",
           "CHEMBL287556": "kojic acid", "CHEMBL232202": "arbutin",
           "CHEMBL165": "resveratrol", "CHEMBL24147": "resorcinol"}
print("\n== anchors present in the gating set? ==")
for cid, name in ANCHORS.items():
    row = g[g.molecule_chembl_id == cid]
    if row.empty:
        print(f"  {name:32s} ABSENT — needs hand-curation")
    else:
        r = row.iloc[0]
        rank = int((g.pAct_median > r.pAct_median).sum()) + 1
        print(f"  {name:32s} pAct {r.pAct_median:5.2f}  rank {rank}/{len(g)}")
print("\n  thiamidol                        ABSENT from ChEMBL — hand-curate (PubChem CID 71543007)")

print("\n== top 15 by potency ==")
print(g.head(15)[["molecule_chembl_id", "pref_name", "pAct_median",
                  "n_measurements", "flag_noisy"]].to_string(index=False))
