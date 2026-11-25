#!/usr/bin/env python3
"""Redock the real ligand (0TR), force the metal question with a tight box,
and discover smina's atom-type vocabulary for a custom Cu-O term."""
import os, re, json, subprocess
import numpy as np

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD, LIG = f"{ROOT}/models", f"{ROOT}/ligands"
TPL = f"{MOD}/2Y9X.pdb"

het = {}
for l in open(TPL):
    if l.startswith("HETATM") and l[21] == "A":
        het.setdefault(l[17:20].strip(), []).append(l)
SKIP = {"CU", "ZN", "HOH", "HO", "NAG", "GOL", "EDO", "SO4", "PO4", "CL", "NA"}
LIGCODE = max((k for k in het if k not in SKIP), key=lambda k: len(het[k]))
print(f"co-crystal ligand: {LIGCODE} ({len(het[LIGCODE])} atoms)")

cry = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in het[LIGCODE]])
cry_o = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                  for l in het[LIGCODE] if l[76:78].strip() == "O"])
cu = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in het["CU"]])[:2]
c = cry.mean(0)
print(f"  crystal centroid {c.round(2)}   min O-Cu {min(np.linalg.norm(o-m) for o in cry_o for m in cu):.2f} A")

print("\n== how are the coppers written in the receptor pdbqt? ==")
for l in open(f"{MOD}/hTYR_receptor.pdbqt"):
    if l.startswith(("ATOM", "HETATM")) and "CU" in l[12:20].upper():
        print("   ", repr(l.rstrip()))

rec = f"{MOD}/2Y9X_rec.pdb"
with open(rec, "w") as f:
    for l in open(TPL):
        if l[21:22] != "A": continue
        if l.startswith("ATOM") or (l.startswith("HETATM") and l[17:20].strip() == "CU"):
            f.write(l)
    f.write("END\n")
open(f"{MOD}/{LIGCODE}.pdb", "w").write("".join(het[LIGCODE]) + "END\n")
os.system(f"obabel {rec} -O {MOD}/2Y9X_rec.pdbqt -xr 2>/dev/null")
os.system(f"obabel {MOD}/{LIGCODE}.pdb -O {MOD}/{LIGCODE}.pdbqt -h 2>/dev/null")

nocu = f"{MOD}/2Y9X_rec_nocu.pdbqt"
with open(nocu, "w") as f:
    for l in open(f"{MOD}/2Y9X_rec.pdbqt"):
        if l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu"): continue
        f.write(l)

def run(rec, lig, cx, cy, cz, size, extra=()):
    out = f"{MOD}/_rd.pdbqt"
    p = subprocess.run(["smina", "-r", rec, "-l", lig,
        "--center_x", f"{cx:.2f}", "--center_y", f"{cy:.2f}", "--center_z", f"{cz:.2f}",
        "--size_x", str(size), "--size_y", str(size), "--size_z", str(size),
        "--exhaustiveness", "16", "--num_modes", "3", "--seed", "42",
        "--cpu", "4", "-o", out, *extra], capture_output=True, text=True)
    m = re.search(r"^\s*1\s+(-?\d+\.\d+)", p.stdout, re.M)
    if not m:
        return None, None, None, (p.stderr or p.stdout)[-200:]
    xyz, o, on = [], [], False
    for l in open(out):
        if l.startswith("MODEL"): on = True
        elif on and l.startswith(("ATOM", "HETATM")):
            v = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            xyz.append(v)
            if l.split()[-1] in ("OA", "O"): o.append(v)
        elif l.startswith("ENDMDL") and xyz: break
    return float(m.group(1)), np.array(xyz), (np.array(o) if o else None), ""

print(f"\n== redocking {LIGCODE} into its own crystal structure (12 A box) ==")
for tag, r in (("Cu present", f"{MOD}/2Y9X_rec.pdbqt"), ("no metal", nocu)):
    sc, xyz, o, err = run(r, f"{MOD}/{LIGCODE}.pdbqt", *c, 12)
    if sc is None: print(f"  {tag:12s} FAILED: {err}"); continue
    disp = np.linalg.norm(xyz.mean(0) - c)
    ocu = min(np.linalg.norm(a-b) for a in o for b in cu) if o is not None else float("nan")
    print(f"  {tag:12s} score {sc:6.2f}   centroid off {disp:5.2f} A   min O-Cu {ocu:5.2f} A")

print("\n== tight box forces the metal question (human model, 10 A) ==")
box = json.load(open(f"{MOD}/box.json"))
for tag, r in (("Cu present", f"{MOD}/hTYR_receptor.pdbqt"),
               ("Cu as Zn",   f"{MOD}/hTYR_receptor_zn.pdbqt"),
               ("no metal",   f"{MOD}/hTYR_receptor_nocu.pdbqt")):
    sc, xyz, o, err = run(r, f"{LIG}/CHEMBL537.pdbqt",
                          box["center_x"], box["center_y"], box["center_z"], 10)
    print(f"  hydroquinone, {tag:11s}: {sc if sc is not None else 'FAILED ' + err}")

print("\n== probing smina's atom-type names ==")
BASE = """-0.035579 gauss(o=0,_w=0.5,_c=8)
-0.005156 gauss(o=3,_w=2,_c=8)
0.840245 repulsion(o=0,_c=8)
-0.035069 hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439 non_dir_h_bond(g=-0.7,_b=0,_c=8)
"""
for name in ("Copper", "Cu", "Zinc", "GenericMetal", "Magnesium", "Iron"):
    f = f"{MOD}/_probe.score"
    open(f, "w").write(BASE + f"-1.0 atom_type_gaussian(t1={name},t2=Oxygen,o=2.2,_w=0.8,_c=8)\n")
    p = subprocess.run(["smina", "-r", f"{MOD}/hTYR_receptor.pdbqt",
        "-l", f"{LIG}/CHEMBL537.pdbqt", "--custom_scoring", f, "--score_only"],
        capture_output=True, text=True)
    ok = "Affinity" in p.stdout
    msg = "" if ok else (p.stderr or p.stdout).strip().splitlines()[:1]
    print(f"  t1={name:14s} {'ACCEPTED' if ok else 'rejected'}  {msg if msg else ''}")
