#!/usr/bin/env python3
"""For the 19 'no evidence' documents: pull metadata and abstracts, then
classify which enzyme they actually used."""
import os, re, time
import pandas as pd, requests
from chembl_webresource_client.new_client import new_client

ROOT = os.path.expanduser("~/hTYR_Platform")
RAW, CUR = f"{ROOT}/data/raw", f"{ROOT}/data/curated"
tri = pd.read_csv(f"{CUR}/human_tyr_triaged.csv")
allact = pd.read_csv(f"{RAW}/tyrosinase_activities.csv")

fung = {"CHEMBL3318", "CHEMBL3988582", "CHEMBL6066154"}
fdocs = set(allact[allact.target_chembl_id.isin(fung)].document_chembl_id.dropna())
prof = tri.groupby("document_chembl_id").bucket.value_counts().unstack(fill_value=0)
for c in ("human_enzyme", "human_cell", "unclear", "suspect_fungal"):
    if c not in prof: prof[c] = 0
docs = [d for d, r in prof.iterrows()
        if r.unclear > 0 and r.human_enzyme == 0 and r.human_cell == 0
        and r.suspect_fungal == 0 and d not in fdocs]
print(f"{len(docs)} documents to check\n")

doc_api = new_client.document
FUNGAL = re.compile(r"mushroom|agaricus|fungal|polyphenol oxidase|bisporus", re.I)
HUMAN  = re.compile(r"human tyrosinase|recombinant human|hTYR|HEMn|MNT-?1|"
                    r"human melanocyte|B16|murine|mouse|melan-a", re.I)
rows = []
for d in docs:
    try:
        m = doc_api.get(d)
    except Exception as e:
        print(f"  {d}: lookup failed ({e})"); continue
    pmid, doi = m.get("pubmed_id"), m.get("doi")
    title = (m.get("title") or "")[:95]
    text = title
    if pmid:
        try:
            r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                             params={"db": "pubmed", "id": pmid, "rettype": "abstract",
                                     "retmode": "text"},
                             headers={"User-Agent": "htyr-benchmark/0.1"}, timeout=45)
            if r.ok: text += " " + r.text
            time.sleep(0.4)
        except Exception:
            pass
    f_hit = sorted(set(w.lower() for w in FUNGAL.findall(text)))
    h_hit = sorted(set(w.lower() for w in HUMAN.findall(text)))
    call = ("FUNGAL" if f_hit and not h_hit else
            "HUMAN"  if h_hit and not f_hit else
            "BOTH"   if f_hit and h_hit else "UNKNOWN")
    n = int(prof.loc[d, "unclear"])
    rows.append({"document_chembl_id": d, "pmid": pmid, "doi": doi, "title": title,
                 "n_unclear": n, "fungal_terms": ";".join(f_hit),
                 "human_terms": ";".join(h_hit), "call": call})
    print(f"  {d} [{call:7s}] n={n:3d}  {title}")
    if f_hit or h_hit:
        print(f"      fungal={f_hit}  human={h_hit}")

out = pd.DataFrame(rows)
out.to_csv(f"{CUR}/unclear_doc_calls.csv", index=False)
print("\n== summary ==")
for c, g in out.groupby("call"):
    print(f"  {c:8s} {len(g):3d} documents, {int(g.n_unclear.sum()):4d} assays")
print(f"\nwrote {CUR}/unclear_doc_calls.csv")
print("Anything called HUMAN still needs the methods section read before use,")
print("and carries weaker provenance than the locked 74 — flag it, don't merge silently.")
