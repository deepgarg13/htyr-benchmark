#!/usr/bin/env python3
"""Fetch the AlphaFold model for human TYR and candidate metal-donor templates,
then report which metals each template actually contains and what coordinates them."""
import os, requests

MOD = os.path.expanduser("~/hTYR_Platform/models")
os.makedirs(MOD, exist_ok=True)

print("== AlphaFold model for P14679 ==")
meta = requests.get("https://alphafold.ebi.ac.uk/api/prediction/P14679", timeout=60).json()[0]
for key, ext in (("pdbUrl", "pdb"), ("cifUrl", "cif")):
    if meta.get(key):
        r = requests.get(meta[key], timeout=120)
        open(f"{MOD}/AF_P14679.{ext}", "wb").write(r.content)
        print(f"  {meta[key].split('/')[-1]}  ({len(r.content)//1024} KB)")
print(f"  model version: {meta.get('latestVersion')}  seq len: {len(meta.get('uniprotSequence',''))}")

TEMPLATES = {
    "2Y9X": "Agaricus bisporus tyrosinase (mushroom) — Cu donor",
    "3NM8": "Bacillus megaterium tyrosinase — Cu donor",
    "5M8P": "Human TYRP1 — closest human fold (check its metal!)",
}
for pdb, note in TEMPLATES.items():
    r = requests.get(f"https://files.rcsb.org/download/{pdb}.pdb", timeout=120)
    if r.status_code != 200:
        print(f"\n{pdb}: download failed ({r.status_code})"); continue
    open(f"{MOD}/{pdb}.pdb", "wb").write(r.content)
    lines = r.text.splitlines()
    metals = [l for l in lines if l.startswith("HETATM")
              and l[17:20].strip() in ("CU", "ZN", "FE", "MN", "NI")]
    print(f"\n{pdb}: {note}")
    print(f"  metal atoms found: {len(metals)}")
    coords = []
    for m in metals:
        el = m[17:20].strip()
        xyz = (float(m[30:38]), float(m[38:46]), float(m[46:54]))
        coords.append((el, xyz))
        print(f"    {el:3s} chain {m[21]} resseq {m[22:26].strip():>5s}")
    # metal-metal distances
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            d = sum((a - b) ** 2 for a, b in zip(coords[i][1], coords[j][1])) ** 0.5
            if d < 6.0:
                print(f"    {coords[i][0]}–{coords[j][0]} distance: {d:.2f} A")
    # coordinating residues within 2.6 A
    prot = [l for l in lines if l.startswith("ATOM")]
    for el, xyz in coords[:2]:
        near = set()
        for a in prot:
            ax = (float(a[30:38]), float(a[38:46]), float(a[46:54]))
            d = sum((p - q) ** 2 for p, q in zip(xyz, ax)) ** 0.5
            if d < 2.6:
                near.add(f"{a[17:20].strip()}{a[22:26].strip()}:{a[12:16].strip()}")
        print(f"    {el} coordinated by: {sorted(near)}")
