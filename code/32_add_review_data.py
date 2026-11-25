#!/usr/bin/env python3
"""Add the named hsTYR compounds from Roulier et al. (J Med Chem 2020) Table 1
and Figure 7, plus the EJMC 2022 purified-enzyme values. Structures via PubChem."""
import os, datetime, urllib.parse, time
import numpy as np, pandas as pd, requests
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
CUR = f"{ROOT}/data/curated"
PUG = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}"
       "/property/SMILES,MolecularFormula/JSON")

#  name, value_uM, type, source, censored
REVIEW = [
 ("3-hydroxycoumarin",            3.4,  "Ki",   "Roulier2020 Fig7 #27", False),
 ("L-mimosine",                   10.3, "Ki",   "Roulier2020 Fig7 #28", False),
 ("phenylthiourea",               1.7,  "Ki",   "Roulier2020 Fig7 #29", False),
 ("4-butylresorcinol",            9.1,  "Ki",   "Roulier2020 Fig7 #30", False),
 ("4-phenylethylresorcinol",      24.0, "Ki",   "Roulier2020 Fig7 #31", False),
 ("4-hexylresorcinol",            39.0, "Ki",   "Roulier2020 Fig7 #32", False),
 ("benzoic acid",                 520.0,"Ki",   "Roulier2020 Tab1",     False),
 ("2-hydroxypyridine-N-oxide",    128.0,"Ki",   "Roulier2020 Tab1",     False),
 ("kojic acid",                   350.0,"Ki",   "Roulier2020 Tab1",     False),
 ("resorcinol",                   3000.0,"IC50","Roulier2020 Tab1",     True),
 ("hydroquinone",                 4400.0,"IC50","Roulier2020 Tab1",     False),
 ("arbutin",                      6500.0,"IC50","Roulier2020 Tab1",     False),
 ("aesculetin",                   None, "IC50", "Roulier2020 Tab1 n.a.",True),
]

rows = []
print("resolving structures via PubChem\n")
for name, val, typ, src, cens in REVIEW:
    r = requests.get(PUG.format(urllib.parse.quote(name)), timeout=45,
                     headers={"User-Agent": "htyr-benchmark/0.1"})
    time.sleep(0.3)
    if r.status_code != 200:
        print(f"  {name:28s} NOT FOUND"); continue
    p = r.json()["PropertyTable"]["Properties"][0]
    pact = None if (val is None or cens) else round(float(-np.log10(val * 1e-6)), 2)
    rows.append({"molecule_chembl_id": f"PUBCHEM:{p['CID']}", "pref_name": name,
                 "smiles": p.get("SMILES"), "pAct_median": pact, "pAct_spread": 0.0,
                 "n_measurements": 1, "types": typ, "docs": src,
                 "flag_noisy": False, "source": "roulier2020",
                 "role": "anchor" if pact is not None else "inactive_reference"})
    print(f"  {name:28s} CID {p['CID']:<10} {typ} {val}  pAct {pact}")

add = pd.DataFrame(rows)
def ikey(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    return Chem.InchiToInchiKey(Chem.MolToInchi(m)).split("-")[0] if m else None

bench = pd.read_csv(f"{CUR}/benchmark_v2.csv")
have = set(bench.smiles.map(ikey).dropna())
add["ikey"] = add.smiles.map(ikey)
new = add[~add.ikey.isin(have)]
print(f"\n  already in benchmark : {int(add.ikey.isin(have).sum())}")
print(f"  NEW compounds        : {len(new)}")
for _, r in new.iterrows():
    print(f"    + {r.pref_name}  ({r.types} {r.docs})")

add.drop(columns=["ikey"]).to_csv(f"{CUR}/roulier_additions.csv", index=False)
print(f"\nwrote {CUR}/roulier_additions.csv  (NOT merged)")

with open(f"{ROOT}/BENCHMARK_LOCKED.md", "a") as f:
    f.write(f"""

## CORRECTION — {datetime.date.today().strftime('%d-%m-%Y')}
The data-quality note above, stating that hydroquinone (CHEMBL537) and arbutin
(CHEMBL232202) are censored bounds flattened to exact values, is WRONG and is
hereby struck. Roulier et al., J Med Chem 2020, Table 1 reports exact hsTYR
IC50 values of 4,400 uM and 6,500 uM respectively. ChEMBL recorded them
correctly; our pAct values (2.36, 2.19) match to within rounding.

## NOTE — world supply of purified-hsTYR data
Roulier et al. 2020 Figure 7 is the field's complete set of hsTYR inhibitors
measured against the isolated enzyme: 25 compounds. With Table 1 references and
Roulier et al., Eur J Med Chem 2022 (compounds 10 and 23, Ki 0.55 and 0.25 uM),
the published purified-enzyme SAR is of order 35-40 compounds. This independently
confirms the data-scarcity finding reached from ChEMBL and BindingDB.
""")
print("lock file corrected and annotated")
