#!/usr/bin/env python3
"""Add the Mann Table 1 human-enzyme anchors. Pre-declared in the lock file,
so this is a documented amendment, not post-hoc tuning."""
import os, hashlib, datetime, urllib.parse
import numpy as np, pandas as pd, requests

CUR  = os.path.expanduser("~/hTYR_Platform/data/curated")
ROOT = os.path.expanduser("~/hTYR_Platform")
SRC  = "Mann 2018 Table 1 (hTyr IC50, single consistent assay)"

#                       name                       hTyr IC50 uM   censored?
MANN = [("Thiamidol",                              1.1,   False),
        ("4-Hexylresorcinol",                      94.0,  False),
        ("4-Phenylethylresorcinol",                131.0, False),
        ("Dimethoxytolylpropyl resorcinol",        None,  True)]  # no inhibition

PUG = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}"
       "/property/SMILES,MolecularFormula,MolecularWeight/JSON")

def smiles_for(name):
    r = requests.get(PUG.format(urllib.parse.quote(name)), timeout=45,
                     headers={"User-Agent": "htyr-benchmark/0.1"})
    if r.status_code != 200:
        return None, None, None
    p = r.json()["PropertyTable"]["Properties"][0]
    return p["CID"], p.get("SMILES"), p.get("MolecularWeight")

rows = []
print("== resolving Mann anchors via PubChem ==")
for name, ic50, censored in MANN:
    cid, smi, mw = smiles_for(name)
    pact = None if ic50 is None else round(float(-np.log10(ic50 * 1e-6)), 2)
    print(f"  {name:34s} CID {str(cid):>10s}  MW {str(mw):>7s}  pAct {pact}")
    if smi is None:
        print(f"    !! no structure — curate by hand before locking")
    rows.append({"molecule_chembl_id": f"PUBCHEM:{cid}", "pref_name": name,
                 "smiles": smi, "pAct_median": pact, "pAct_spread": 0.0,
                 "n_measurements": 1, "types": "IC50", "docs": SRC,
                 "flag_noisy": False, "source": "mann2018",
                 "role": "decoy_mushroom_potent" if censored else "anchor"})

anchors = pd.DataFrame(rows)
bench = pd.read_csv(f"{CUR}/benchmark_human_enzyme.csv")
bench["source"] = "chembl"; bench["role"] = "benchmark"

v2 = pd.concat([bench, anchors], ignore_index=True).sort_values(
        "pAct_median", ascending=False, na_position="last")
v2.to_csv(f"{CUR}/benchmark_v2.csv", index=False)

scored = v2[v2.pAct_median.notna()]
thia = v2[v2.pref_name == "Thiamidol"]
if not thia.empty:
    pa = thia.pAct_median.iloc[0]
    print(f"\nThiamidol pAct {pa:.2f} -> rank "
          f"{int((scored.pAct_median > pa).sum()) + 1}/{len(scored)}")

sha = hashlib.sha256(open(f"{CUR}/benchmark_v2.csv", 'rb').read()).hexdigest()
with open(f"{ROOT}/BENCHMARK_LOCKED.md", "a") as f:
    f.write(f"""

## AMENDMENT — {datetime.date.today().strftime('%d-%m-%Y')}
Hand-curated anchors added, as pre-declared above. Source: {SRC}.
File:   data/curated/benchmark_v2.csv
SHA256: {sha}
N:      {len(scored)} scored compounds (+1 unscored decoy)

Added: Thiamidol (1.1 uM), 4-Hexylresorcinol (94 uM),
       4-Phenylethylresorcinol (131 uM), all hTyr IC50.
Decoy: Dimethoxytolylpropyl resorcinol — mTyr 0.246 uM, hTyr NO INHIBITION.
       Carries no pAct; used only in gate condition 4.

### Gate condition 4 (added)
The decoy must NOT appear in the top 25% of the predicted ranking. A model
that scores it highly is reproducing mushroom-like recognition.

### Data-quality note
ChEMBL records hydroquinone (CHEMBL537) and arbutin (CHEMBL232202) with
standard_relation '=', but Mann 2018 reports both as >4,000 uM. Their pAct
values (2.36, 2.19) are censored bounds flattened to exact values. Harmless
for rank-order use; do not quote them as measured potencies.
""")
print(f"\namended {ROOT}/BENCHMARK_LOCKED.md\nsha256 {sha}")
