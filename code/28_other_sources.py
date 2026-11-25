#!/usr/bin/env python3
"""Is there human-TYR potency data outside ChEMBL? BindingDB and PubChem."""
import os, json, time
import pandas as pd, requests

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
UA = {"User-Agent": "htyr-benchmark/0.1"}
ACC = "P14679"

print("== BindingDB ==")
for url in (f"https://bindingdb.org/rest/getLigandsByUniprots?uniprot={ACC}"
            f"&cutoff=10000000&response=application/json",
            f"https://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot"
            f"?uniprot={ACC}"):
    try:
        r = requests.get(url, timeout=90, headers=UA)
        print(f"  {url.split('?')[0].split('/')[-1]}: HTTP {r.status_code}, "
              f"{len(r.content)} bytes, {r.headers.get('content-type','?')}")
        if r.ok and len(r.content) > 200:
            head = r.text[:400].replace("\n", " ")
            print(f"    {head}")
            open(f"{CUR}/_bindingdb_raw.txt", "w").write(r.text)
            print(f"    saved -> {CUR}/_bindingdb_raw.txt")
            break
    except Exception as e:
        print(f"  failed: {e}")

print("\n== PubChem BioAssay ==")
try:
    r = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/"
                     f"accession/{ACC}/aids/JSON", timeout=60, headers=UA)
    print(f"  aids lookup: HTTP {r.status_code}")
    if r.ok:
        aids = r.json()["IdentifierList"]["AID"]
        print(f"  {len(aids)} assays reference {ACC}")
        rows = []
        for aid in aids[:40]:
            s = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"
                             f"{aid}/summary/JSON", timeout=45, headers=UA)
            time.sleep(0.25)
            if not s.ok: continue
            d = s.json()["AssaySummaries"]["AssaySummary"][0]
            name = (d.get("Name") or "")[:80]
            n_act = d.get("ActiveCount", 0); n_all = d.get("CIDCountAll", 0)
            rows.append({"aid": aid, "name": name, "active": n_act, "tested": n_all})
            print(f"    AID {aid:>8}  tested={n_all:<7} active={n_act:<6} {name}")
        if rows:
            pd.DataFrame(rows).to_csv(f"{CUR}/_pubchem_assays.csv", index=False)
            print(f"\n  saved -> {CUR}/_pubchem_assays.csv")
    else:
        print(f"  {r.text[:200]}")
except Exception as e:
    print(f"  failed: {e}")

print("\nLook for assays with dose-response (IC50) rather than single-point screens,")
print("and check each one's description names the HUMAN enzyme before using it.")
