#!/usr/bin/env python3
"""Feasibility gate: how much potency data exists for HUMAN tyrosinase
(UniProt P14679) versus the mushroom enzyme the field actually screens?"""
import os
import pandas as pd
from chembl_webresource_client.new_client import new_client

OUT = os.path.expanduser("~/hTYR_Platform/data/raw")
os.makedirs(OUT, exist_ok=True)

HUMAN_ACC = "P14679"          # TYRO_HUMAN
POTENCY = {"IC50", "Ki", "EC50", "Kd"}
ANCHORS = ["isobutylamido thiazolyl resorcinol", "4-butylresorcinol",
           "rucinol", "cysteamine", "hydroquinone", "kojic acid",
           "arbutin", "tranexamic acid", "resveratrol", "azelaic acid"]

target, activity, molecule = new_client.target, new_client.activity, new_client.molecule

ACT_FIELDS = ["molecule_chembl_id", "standard_type", "standard_relation",
              "standard_value", "standard_units", "pchembl_value",
              "assay_description", "document_chembl_id"]

def pull(tid):
    df = pd.DataFrame(list(activity.filter(target_chembl_id=tid).only(ACT_FIELDS)))
    if not df.empty:
        df["target_chembl_id"] = tid
    return df

print("== human tyrosinase targets ==")
human = list(target.filter(target_components__accession=HUMAN_ACC))
for t in human:
    print(f"  {t['target_chembl_id']}  {t['pref_name']}  [{t['target_type']}]")

print("\n== mushroom / fungal tyrosinase targets ==")
fungal = [t for t in target.filter(organism__icontains="Agaricus")]
for t in fungal:
    print(f"  {t['target_chembl_id']}  {t['pref_name']}  ({t['organism']})")

frames, counts = [], {}
for label, tl in (("human", human), ("fungal", fungal)):
    for t in tl:
        tid = t["target_chembl_id"]
        df = pull(tid)
        if df.empty:
            print(f"\n{label} {tid}: no activities"); continue
        df["group"] = label
        pot = df[df.standard_type.isin(POTENCY)]
        counts.setdefault(label, [0, 0])
        counts[label][0] += len(pot)
        counts[label][1] += pot.molecule_chembl_id.nunique()
        print(f"\n{label} {tid}: {len(df)} activities, {len(pot)} potency, "
              f"{pot.molecule_chembl_id.nunique()} distinct compounds")
        frames.append(df)

if frames:
    allact = pd.concat(frames, ignore_index=True)
    allact.to_csv(f"{OUT}/tyrosinase_activities.csv", index=False)
    print(f"\nwrote {len(allact)} rows -> {OUT}/tyrosinase_activities.csv")

print("\n== THE FEASIBILITY NUMBER ==")
for label in ("human", "fungal"):
    n = counts.get(label, [0, 0])
    print(f"  {label:7s}: {n[0]:6d} potency measurements, {n[1]:5d} distinct compounds")

print("\n== resolving anchor compounds ==")
rows = []
for name in ANCHORS:
    hits = list(molecule.filter(molecule_synonyms__molecule_synonym__icontains=name)
                        .only(["molecule_chembl_id", "pref_name"])[:1])
    if not hits:
        hits = list(molecule.search(name).only(["molecule_chembl_id", "pref_name"])[:1])
    if not hits:
        print(f"  {name:38s} -> NOT FOUND (curate by hand)"); continue
    h = hits[0]
    full = molecule.get(h["molecule_chembl_id"])
    smi = (full.get("molecule_structures") or {}).get("canonical_smiles")
    rows.append({"query": name, "chembl_id": h["molecule_chembl_id"],
                 "pref_name": h.get("pref_name"), "smiles": smi})
    print(f"  {name:38s} -> {h['molecule_chembl_id']}  {h.get('pref_name')}")

pd.DataFrame(rows).to_csv(f"{OUT}/anchors_resolved.csv", index=False)
print("\nDone. Review both CSVs before locking anything.")
