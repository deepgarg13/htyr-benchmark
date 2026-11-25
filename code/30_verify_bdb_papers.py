#!/usr/bin/env python3
"""Which BindingDB source papers actually assayed HUMAN tyrosinase?"""
import os, re, time
import pandas as pd, requests

CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
new = pd.read_csv(f"{CUR}/bindingdb_new.csv")
new["pmid"] = pd.to_numeric(new.pmid, errors="coerce")
new["internal"] = new.doi.astype(str).str.startswith("10.7270")

print(f"{len(new)} new compounds")
print(f"  BindingDB-internal (no paper): {int(new.internal.sum())}")
print(f"  from journal papers          : {int((~new.internal).sum())}")

papers = (new[~new.internal & new.pmid.notna()]
          .assign(pmid=lambda d: d.pmid.astype("int64").astype(str))
          .groupby(["pmid", "doi"]).size().sort_values(ascending=False))
print(f"  checkable papers             : {len(papers)}\n")

FUNGAL = re.compile(r"mushroom|agaricus|bisporus|fungal|polyphenol oxidase", re.I)
HENZ   = re.compile(r"human tyrosinase|recombinant human tyrosin|hTYR\b|human TYR\b", re.I)
HCELL  = re.compile(r"MNT-?1|HEMn|HMVII|SK-?MEL|human melanocyte|human melanoma|"
                    r"human epidermal melanocyte", re.I)
MURINE = re.compile(r"\bB16\b|melan-a|murine|\bmouse\b", re.I)

rows = []
for (pmid, doi), n in papers.items():
    text = ""
    try:
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                         params={"db": "pubmed", "id": pmid, "rettype": "abstract",
                                 "retmode": "text"},
                         headers={"User-Agent": "htyr-benchmark/0.1"}, timeout=45)
        if r.ok: text = r.text
        time.sleep(0.4)
    except Exception as e:
        print(f"  PMID {pmid}: fetch failed ({e})")
    f, he, hc, mu = (bool(p.search(text)) for p in (FUNGAL, HENZ, HCELL, MURINE))
    call = ("FUNGAL" if f and not (he or hc) else
            "HUMAN_ENZYME" if he and not f else
            "HUMAN_CELL" if hc and not (f or he) else
            "MIXED" if (f and (he or hc)) else
            "MURINE" if mu and not (f or he or hc) else "UNKNOWN")
    title = next((l.strip() for l in text.splitlines()
                  if len(l.strip()) > 40 and not l.startswith(("1.", "Author"))), "")[:85]
    rows.append({"pmid": pmid, "doi": doi, "n_compounds": int(n), "call": call,
                 "title": title})
    print(f"  PMID {pmid:>9} [{call:12s}] {n:3d} cmpds  {title}")

res = pd.DataFrame(rows)
if res.empty:
    raise SystemExit("no papers checked - inspect the pmid column of bindingdb_new.csv")
res.to_csv(f"{CUR}/bdb_paper_calls.csv", index=False)
print("\n== recoverable by provenance class ==")
for c, g in res.groupby("call"):
    print(f"  {c:13s} {len(g):3d} papers, {int(g.n_compounds.sum()):4d} compounds")

usable = res[res.call == "HUMAN_ENZYME"]
print(f"\nDEFENSIBLE human-enzyme compounds available: "
      f"{int(usable.n_compounds.sum())} from {len(usable)} papers")
print(f"wrote {CUR}/bdb_paper_calls.csv")
