#!/usr/bin/env python3
"""Are the two key papers open access? Fetch what's legally available."""
import os, requests, time

OUT = os.path.expanduser("~/hTYR_Platform/papers")
os.makedirs(OUT, exist_ok=True)
UA = {"User-Agent": "htyr-benchmark/0.1"}
WANT = {"36462443": "resorcinol hemiindigoids (22 cmpds)",
        "32787103": "Roulier - genuine human TYR inhibitors (11 cmpds)",
        "16392817": "aurones (3 cmpds, classed HUMAN_ENZYME)"}

for pmid, note in WANT.items():
    print(f"\nPMID {pmid} - {note}")
    r = requests.get("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
                     params={"ids": pmid, "format": "json"}, headers=UA, timeout=45)
    rec = (r.json().get("records") or [{}])[0] if r.ok else {}
    pmcid, doi = rec.get("pmcid"), rec.get("doi")
    print(f"  doi={doi}  pmcid={pmcid}")
    if not pmcid:
        print("  not in PMC - paywalled. Read via your institution or request a copy.")
        continue
    f = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                     params={"db": "pmc", "id": pmcid.replace("PMC", ""),
                             "rettype": "full", "retmode": "xml"},
                     headers=UA, timeout=90)
    if f.ok and len(f.content) > 5000:
        p = f"{OUT}/{pmid}_{pmcid}.xml"
        open(p, "wb").write(f.content)
        print(f"  saved {len(f.content)//1024} KB -> {p}")
    else:
        print(f"  PMC returned {f.status_code}, {len(f.content)} bytes "
              f"(likely not in the open-access subset)")
    time.sleep(0.5)

print("\nFor anything paywalled: the data you need is the human-TYR IC50 column")
print("of the SAR table. Paste it here like you did with Mann Table 1 and I'll")
print("curate it into the benchmark with provenance.")
