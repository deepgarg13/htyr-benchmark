#!/usr/bin/env python3
"""Redocking control in the native structure, and whether retyping the metal
makes smina see it. Attempt 3."""
import os, re, subprocess
import numpy as np

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD, LIG = f"{ROOT}/models", f"{ROOT}/ligands"
TPL = f"{MOD}/2Y9X.pdb"

print("== what terms does this smina build offer? ==")
p = subprocess.run(["smina", "--print_terms"], capture_output=True, text=True)
for l in (p.stdout + p.stderr).splitlines():
    if re.search(r"metal|zinc|atom_type|ad4", l, re.I):
        print("   ", l.strip())

print("\n== het groups in 2Y9X chain A ==")
het = {}
for l in open(TPL):
    if l.startswith("HETATM") and l[21] == "A":
        het.setdefault(l[17:20].strip(), []).append(l)
for k, v in het.items():
    print(f"  {k:5s} {len(v):3d} atoms")

LIGCODE = next((k for k in het if k not in ("CU", "HOH", "ZN", "NAG", "GOL", "EDO")), None)
print(f"  -> treating {LIGCODE} as the co-crystal ligand")

cry = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                for l in het[LIGCODE] if l[76:78].strip() != "H"])
cry_o = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                  for l in het[LIGCODE] if l[76:78].strip() == "O"])
cu = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in het["CU"]])[:2]
print(f"  crystal pose: centroid {cry.mean(0).round(2)}")
if len(cry_o):
    print(f"  crystal min O-Cu distance: {min(np.linalg.norm(o-c) for o in cry_o for c in cu):.2f} A")

# receptor from 2Y9X chain A: protein + Cu, ligand and waters removed
rec_pdb = f"{MOD}/2Y9X_rec.pdb"
with open(rec_pdb, "w") as f:
    for l in open(TPL):
        if l[21:22] != "A": continue
        if l.startswith("ATOM") or (l.startswith("HETATM") and l[17:20].strip() == "CU"):
            f.write(l)
    f.write("END\n")
lig_pdb = f"{MOD}/{LIGCODE}.pdb"
open(lig_pdb, "w").write("".join(het[LIGCODE]) + "END\n")
os.system(f"obabel {rec_pdb} -O {MOD}/2Y9X_rec.pdbqt -xr 2>/dev/null")
os.system(f"obabel {lig_pdb} -O {MOD}/{LIGCODE}.pdbqt -h 2>/dev/null")

def retype(src, dst, new="Zn"):
    with open(dst, "w") as f:
        for l in open(src):
            if l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu"):
                l = l.rstrip()[:-2] + new + "\n"
            f.write(l)

retype(f"{MOD}/2Y9X_rec.pdbqt", f"{MOD}/2Y9X_rec_zn.pdbqt")
retype(f"{MOD}/hTYR_receptor.pdbqt", f"{MOD}/hTYR_receptor_zn.pdbqt")
nocu = f"{MOD}/2Y9X_rec_nocu.pdbqt"
with open(nocu, "w") as f:
    for l in open(f"{MOD}/2Y9X_rec.pdbqt"):
        if l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu"): continue
        f.write(l)

c = cry.mean(0)
def dock(rec, out):
    subprocess.run(["smina", "-r", rec, "-l", f"{MOD}/{LIGCODE}.pdbqt",
        "--center_x", f"{c[0]:.2f}", "--center_y", f"{c[1]:.2f}",
        "--center_z", f"{c[2]:.2f}", "--size_x", "18", "--size_y", "18",
        "--size_z", "18", "--exhaustiveness", "16", "--num_modes", "5",
        "--seed", "42", "--cpu", "4", "-o", out], capture_output=True, text=True)
    xyz, o, sc, on = [], [], None, False
    for l in open(out):
        if l.startswith("MODEL"): on = True
        elif l.startswith("REMARK"):
            m = re.search(r"(?:minimizedAffinity|VINA RESULT:)\s+(-?\d+\.\d+)", l)
            if m and sc is None: sc = float(m.group(1))
        elif on and l.startswith(("ATOM", "HETATM")):
            v = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            xyz.append(v)
            if l.split()[-1] in ("OA", "O"): o.append(v)
        elif l.startswith("ENDMDL") and xyz: break
    return sc, np.array(xyz), (np.array(o) if o else None)

print(f"\n== redocking {LIGCODE} into its own structure ==")
for tag, rec in (("Cu as written", f"{MOD}/2Y9X_rec.pdbqt"),
                 ("Cu retyped Zn", f"{MOD}/2Y9X_rec_zn.pdbqt"),
                 ("no metal",      nocu)):
    sc, xyz, o = dock(rec, f"{MOD}/_redock_{tag.split()[-1]}.pdbqt")
    disp = np.linalg.norm(xyz.mean(0) - c) if len(xyz) else float("nan")
    ocu = min(np.linalg.norm(a - b) for a in o for b in cu) if o is not None else float("nan")
    print(f"  {tag:15s} score {str(sc):>7s}  centroid off by {disp:5.2f} A"
          f"  min O-Cu {ocu:5.2f} A")
print("  (crystal is the reference: centroid offset ~0 and O-Cu ~2-3 A means the setup works)")

print("\n== does retyping make smina see the metal in the human model? ==")
for tag, rec in (("Cu as written", f"{MOD}/hTYR_receptor.pdbqt"),
                 ("Cu retyped Zn", f"{MOD}/hTYR_receptor_zn.pdbqt"),
                 ("no metal",      f"{MOD}/hTYR_receptor_nocu.pdbqt")):
    import json; box = json.load(open(f"{MOD}/box.json"))
    p = subprocess.run(["smina", "-r", rec, "-l", f"{LIG}/CHEMBL537.pdbqt",
        "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]), "--size_x", "18", "--size_y", "18",
        "--size_z", "18", "--exhaustiveness", "8", "--num_modes", "1",
        "--seed", "42", "--cpu", "4"], capture_output=True, text=True)
    m = re.search(r"^\s*1\s+(-?\d+\.\d+)", p.stdout, re.M)
    print(f"  hydroquinone, {tag:15s}: {m.group(1) if m else 'no score'}")
