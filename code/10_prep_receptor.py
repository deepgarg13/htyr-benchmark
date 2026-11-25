#!/usr/bin/env python3
"""Trim hTYR to its lumenal catalytic domain, protonate, re-seat the coppers,
and define the docking box from the Cu midpoint."""
import os, json, requests
import numpy as np

MOD = os.path.expanduser("~/hTYR_Platform/models")
SRC = f"{MOD}/hTYR_Cu.pdb"

j = requests.get("https://rest.uniprot.org/uniprotkb/P14679.json", timeout=60).json()
print("== topology ==")
bounds = {}
for f in j.get("features", []):
    t = f.get("type")
    if t in ("Signal", "Transmembrane", "Chain", "Topological domain", "Domain"):
        a = f["location"]["start"]["value"]; b = f["location"]["end"]["value"]
        print(f"  {t:20s} {a:4d}-{b:4d}  {f.get('description','')}")
        bounds.setdefault(t, []).append((a, b, f.get("description", "")))

sig = bounds.get("Signal", [(0, 18, "")])[0][1]
tm  = bounds.get("Transmembrane", [])
start, end = sig + 1, (tm[0][0] - 1 if tm else 476)
print(f"\n  keeping lumenal domain: {start}-{end}")

cu_lines, prot_lines, plddt = [], [], []
for l in open(SRC):
    if l.startswith("HETATM") and l[17:20].strip() == "CU":
        cu_lines.append(l)
    elif l.startswith("ATOM"):
        r = int(l[22:26])
        if start <= r <= end:
            prot_lines.append(l)
            if l[12:16].strip() == "CA":
                plddt.append(float(l[60:66]))

plddt = np.array(plddt)
print(f"  residues kept: {len(plddt)}")
print(f"  pLDDT  mean {plddt.mean():.1f}  min {plddt.min():.1f}  "
      f"below 70: {(plddt < 70).sum()} residues")

tmp = f"{MOD}/_trimmed_protein.pdb"
open(tmp, "w").write("".join(prot_lines) + "END\n")

print("\n== protonating (pdbfixer, pH 6.8) ==")
from pdbfixer import PDBFixer
from openmm.app import PDBFile
fx = PDBFixer(filename=tmp)
fx.findMissingResidues(); fx.missingResidues = {}      # do not build loops
fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
fx.findMissingAtoms(); fx.addMissingAtoms()
fx.addMissingHydrogens(6.8)
fixed = f"{MOD}/_fixed_protein.pdb"
PDBFile.writeFile(fx.topology, fx.positions, open(fixed, "w"))
print(f"  atoms after H addition: {fx.topology.getNumAtoms()}")

rec = f"{MOD}/hTYR_receptor.pdb"
with open(rec, "w") as f:
    ser = 0
    for l in open(fixed):
        if l.startswith(("ATOM", "HETATM")):
            ser = int(l[6:11]); f.write(l)
    for l in cu_lines:
        ser += 1
        f.write(l[:6] + f"{ser:5d}" + l[11:])
    f.write("END\n")
print(f"  wrote {rec}")

cu = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in cu_lines])
c = cu.mean(0)
print(f"\n== docking box ==")
print(f"  Cu midpoint: {c[0]:.2f} {c[1]:.2f} {c[2]:.2f}   Cu-Cu {np.linalg.norm(cu[0]-cu[1]):.2f} A")
box = dict(center_x=round(float(c[0]), 2), center_y=round(float(c[1]), 2),
           center_z=round(float(c[2]), 2), size_x=24, size_y=24, size_z=24)
json.dump(box, open(f"{MOD}/box.json", "w"), indent=2)
for k, v in box.items():
    print(f"  {k} = {v}")

print("\n== converting to pdbqt and checking Cu typing ==")
os.system(f"obabel {rec} -O {MOD}/hTYR_receptor.pdbqt -xr 2>/dev/null")
cu_q = [l for l in open(f"{MOD}/hTYR_receptor.pdbqt") if " CU" in l[:20].upper()]
print(f"  Cu atoms surviving conversion: {len(cu_q)}")
for l in cu_q:
    print(f"    {l.rstrip()[-30:]}")
if len(cu_q) < 2:
    print("  !! coppers lost in pdbqt conversion — dock against the PDB instead,")
    print("     smina reads it directly and keeps the metals")
