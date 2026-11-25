#!/usr/bin/env python3
"""Triage human-TYR activities by what was actually assayed, and
locate thiamidol, which ChEMBL does not appear to carry."""
import os, re, json, urllib.parse
import requests
import pandas as pd

RAW = os.path.expanduser("~/hTYR_Platform/data/raw")
CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
os.makedirs(CUR, exist_ok=True)

df = pd.read_csv(f"{RAW}/tyrosinase_activities.csv")
h = df[(df.target_chembl_id == "CHEMBL1973")
       & df.standard_type.isin(["IC50", "Ki", "EC50", "Kd"])].copy()
h["assay_description"] = h.assay_description.fillna("")

PATTERNS = [
    ("suspect_fungal", r"mushroom|agaricus|fungal|polyphenol oxidase"),
    ("murine_cell",    r"\bB16|melan-a|murine|mouse"),
    ("human_cell",     r"MNT-1|MNT1|melanocyte|SK-MEL|human.*melanoma|HEMn|normal human"),
    ("human_enzyme",   r"recombinant human|purified human|human tyrosinase|hTYR"),
]

def classify(desc):
    d = desc.lower()
    for label, pat in PATTERNS:
        if re.search(pat, d, re.I):
            return label
    return "unclear"

h["bucket"] = h.assay_description.map(classify)

print("== what was actually assayed (human-annotated target) ==")
for bucket, grp in h.groupby("bucket"):
    print(f"  {bucket:15s} {len(grp):5d} measurements  "
          f"{grp.molecule_chembl_id.nunique():5d} compounds")

usable = h[h.bucket.isin(["human_enzyme", "human_cell"])]
print(f"\n  USABLE (human enzyme or human cell): "
      f"{len(usable)} measurements, {usable.molecule_chembl_id.nunique()} compounds")
print(f"  pchembl present on {usable.pchembl_value.notna().sum()} of them")

h.to_csv(f"{CUR}/human_tyr_triaged.csv", index=False)
print(f"\nwrote {CUR}/human_tyr_triaged.csv")

print("\n== sample assay descriptions per bucket (eyeball these) ==")
for bucket, grp in h.groupby("bucket"):
    print(f"\n--- {bucket} ---")
    for d in grp.assay_description.drop_duplicates().head(4):
        print("   ", d[:150])

print("\n== hunting thiamidol in PubChem ==")
NAMES = ["Isobutylamido thiazolyl resorcinol", "Thiamidol",
         "Isobutylamido-thiazolyl-resorcinol", "Thiamidol 630"]
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/SMILES,MolecularFormula,MolecularWeight,IUPACName/JSON"
for name in NAMES:
    url = PUG.format(urllib.parse.quote(name))
    r = requests.get(url, timeout=45, headers={"User-Agent": "htyr-benchmark/0.1"})
    if r.status_code != 200:
        print(f"  {name:38s} -> not found ({r.status_code})"); continue
    for p in r.json()["PropertyTable"]["Properties"]:
        print(f"  {name:38s} -> CID {p['CID']}  {p.get('MolecularFormula')}  "
              f"MW {p.get('MolecularWeight')}")
        print(f"      {p.get('SMILES')}")
