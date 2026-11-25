#!/usr/bin/env python3
"""Amend the gate; then ask whether the active site is open and whether
smina sees the coppers at all."""
import os, re, datetime, subprocess, json
import numpy as np

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD, LIG = f"{ROOT}/models", f"{ROOT}/ligands"

with open(f"{ROOT}/BENCHMARK_LOCKED.md", "a") as f:
    f.write(f"""

## AMENDMENT — {datetime.date.today().strftime('%d-%m-%Y')} (gate STRENGTHENED)
Measured on dev: pAct ~ heavy-atom count rho = +0.431 (p=1.8e-3, n=50).
A null model ranking by heavy-atom count alone scores rho 0.431 / EF@20% 1.76,
i.e. it PASSES original condition 1. Condition 1 was therefore too permissive.

Condition 1 is REPLACED (strengthened, never weakened):
  1'. Size-corrected Spearman >= 0.30 — rho between pAct and the residual of the
      predicted score regressed on heavy-atom count. Raw rho is still reported.
Conditions 2, 3 and 4 are unchanged. The null model fails 2 (1.76 < 2.0).
Rationale recorded before any model was evaluated against the new condition.
""")
print("gate amended: condition 1 -> size-corrected rho")

def heavy_atoms(path, chain=None):
    out = []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")): continue
        if l[76:78].strip() == "H" or l[12:16].strip().startswith("H"): continue
        if chain and l[21] != chain: continue
        out.append((l[17:20].strip(), l[21], l[22:26].strip(), l[12:16].strip(),
                    np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out

def cu_pair(ats):
    c = [a[4] for a in ats if a[0] == "CU"]
    return np.array(c[:2])

print("\n== is there a cavity next to the coppers? ==")
for label, path, chain in (("AlphaFold hTYR", f"{MOD}/hTYR_receptor.pdb", None),
                           ("2Y9X mushroom", f"{MOD}/2Y9X.pdb", "A")):
    ats = heavy_atoms(path, chain)
    cu = cu_pair(ats)
    if len(cu) < 2:
        print(f"  {label}: coppers not found"); continue
    mid = cu.mean(0)
    prot = [a for a in ats if a[0] not in ("CU", "ZN", "HOH")]
    print(f"\n  {label}  (Cu-Cu {np.linalg.norm(cu[0]-cu[1]):.2f} A)")
    for r in (4, 5, 6, 7, 8):
        n = sum(1 for a in prot if np.linalg.norm(a[4] - mid) < r)
        print(f"    heavy atoms within {r} A of Cu midpoint: {n:4d}")
    d = sorted((np.linalg.norm(a[4] - mid), f"{a[0]}{a[2]}:{a[3]}") for a in prot)
    print(f"    nearest non-metal atoms: " + ", ".join(f"{nm}({dd:.2f})" for dd, nm in d[:6]))

print("\n== does smina actually see the coppers? ==")
nocu = f"{MOD}/hTYR_receptor_nocu.pdbqt"
with open(nocu, "w") as f:
    for l in open(f"{MOD}/hTYR_receptor.pdbqt"):
        if l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu"): continue
        f.write(l)
box = json.load(open(f"{MOD}/box.json"))
probe = f"{LIG}/CHEMBL537.pdbqt"          # hydroquinone
for tag, rec in (("with Cu", f"{MOD}/hTYR_receptor.pdbqt"), ("without Cu", nocu)):
    p = subprocess.run(["smina", "-r", rec, "-l", probe,
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]), "--size_x", "18", "--size_y", "18",
        "--size_z", "18", "--exhaustiveness", "8", "--num_modes", "1",
        "--seed", "42", "--cpu", "4"], capture_output=True, text=True)
    m = re.search(r"^\s*1\s+(-?\d+\.\d+)", p.stdout, re.M)
    print(f"  hydroquinone {tag:12s}: {m.group(1) if m else 'no score'}")
print("  (identical scores => smina is ignoring the metals entirely)")
