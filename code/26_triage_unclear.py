#!/usr/bin/env python3
"""Recover compounds from the 'unclear' bucket by classifying their source
documents, not their assay descriptions."""
import os
import pandas as pd

ROOT = os.path.expanduser("~/hTYR_Platform")
RAW, CUR = f"{ROOT}/data/raw", f"{ROOT}/data/curated"

tri = pd.read_csv(f"{CUR}/human_tyr_triaged.csv")
allact = pd.read_csv(f"{RAW}/tyrosinase_activities.csv")

fungal_targets = {"CHEMBL3318", "CHEMBL3988582", "CHEMBL6066154"}
fungal_docs = set(allact[allact.target_chembl_id.isin(fungal_targets)]
                  .document_chembl_id.dropna())
print(f"documents appearing under a fungal target: {len(fungal_docs)}")

prof = (tri.groupby("document_chembl_id").bucket
          .value_counts().unstack(fill_value=0))
for col in ("human_enzyme", "human_cell", "unclear", "suspect_fungal", "murine_cell"):
    if col not in prof: prof[col] = 0
prof["also_fungal"] = prof.index.isin(fungal_docs)

def classify(r):
    if r.suspect_fungal > 0 or r.also_fungal:       return "ambiguous_or_fungal"
    if r.human_enzyme > 0:                          return "human_enzyme_doc"
    if r.human_cell > 0:                            return "human_cell_doc"
    return "no_evidence"
prof["doc_class"] = prof.apply(classify, axis=1)

print("\n== documents carrying 'unclear' assays ==")
u = prof[prof.unclear > 0]
for cls, grp in u.groupby("doc_class"):
    n_assay = int(grp.unclear.sum())
    print(f"  {cls:22s} {len(grp):3d} documents, {n_assay:4d} unclear assays")

recover = tri[(tri.bucket == "unclear") &
              tri.document_chembl_id.isin(
                  prof[prof.doc_class == "human_enzyme_doc"].index)]
print(f"\nRECOVERABLE as human enzyme: {len(recover)} measurements, "
      f"{recover.molecule_chembl_id.nunique()} compounds")

print("\n== sample descriptions from recoverable documents (sanity-check these) ==")
for doc, grp in list(recover.groupby("document_chembl_id"))[:6]:
    sib = tri[(tri.document_chembl_id == doc) & (tri.bucket == "human_enzyme")]
    print(f"\n  {doc}  ({len(grp)} unclear, {len(sib)} explicitly human)")
    print(f"    explicit : {sib.assay_description.iloc[0][:110] if len(sib) else '-'}")
    print(f"    unclear  : {grp.assay_description.iloc[0][:110]}")

recover.to_csv(f"{CUR}/recoverable_unclear.csv", index=False)
print(f"\nwrote {CUR}/recoverable_unclear.csv")
print("\nNOT merged into the benchmark. Read the samples above first —")
print("if the pairing looks sound, these become an EXTERNAL test set:")
print("never split, never trained on, independent of dev and holdout.")
