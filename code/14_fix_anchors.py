#!/usr/bin/env python3
"""Recover the two anchors PubChem missed, trying synonyms."""
import os, urllib.parse, requests
import pandas as pd

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
PUG = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}"
       "/property/SMILES,MolecularFormula,MolecularWeight/JSON")

TRY = {
 "4-Phenylethylresorcinol": ["4-(1-phenylethyl)benzene-1,3-diol", "Symwhite 377",
                             "phenylethyl resorcinol", "4-(1-phenylethyl)-1,3-benzenediol"],
 "Dimethoxytolylpropyl resorcinol": ["dimethoxytolyl propylresorcinol",
                                     "4-(3-(2,4-dimethoxy-3-methylphenyl)propyl)benzene-1,3-diol",
                                     "dimethoxytolylpropylresorcinol"],
}
found = {}
for canon, names in TRY.items():
    print(f"\n{canon}")
    for n in names:
        r = requests.get(PUG.format(urllib.parse.quote(n)), timeout=45,
                         headers={"User-Agent": "htyr-benchmark/0.1"})
        if r.status_code != 200:
            print(f"  '{n}' -> {r.status_code}"); continue
        p = r.json()["PropertyTable"]["Properties"][0]
        print(f"  '{n}' -> CID {p['CID']}  {p.get('MolecularFormula')}  MW {p.get('MolecularWeight')}")
        print(f"      {p.get('SMILES')}")
        found[canon] = (p["CID"], p.get("SMILES"))
        break

for f in ("benchmark_v2.csv", "dev.csv"):
    path = f"{CUR}/{f}"
    df = pd.read_csv(path)
    for canon, (cid, smi) in found.items():
        m = df.pref_name == canon
        if m.any():
            df.loc[m, "smiles"] = smi
            df.loc[m, "molecule_chembl_id"] = f"PUBCHEM:{cid}"
    df.to_csv(path, index=False)
    print(f"\nupdated {f}")
if len(found) < 2:
    print("\n!! still missing — draw the SMILES by hand from the INCI name before attempt 3")
