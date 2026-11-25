#!/usr/bin/env python3
"""Is the copper in smina's receptor model at all? Errors are shown, not hidden."""
import os, re, subprocess
import numpy as np

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD, LIG = f"{ROOT}/models", f"{ROOT}/ligands"

cu_src = [l for l in open(f"{MOD}/hTYR_receptor.pdbqt")
          if l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu")]
cu = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in cu_src])
mid = cu.mean(0)
base = [l for l in open(f"{MOD}/hTYR_receptor.pdbqt")
        if not (l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu"))]

def variant(path, record, typ):
    with open(path, "w") as f:
        f.writelines(base)
        for l in cu_src:
            f.write(record + l[6:].rstrip()[:-2] + typ + "\n")

V = {"obabel_as_is": f"{MOD}/hTYR_receptor.pdbqt",
     "hetatm_Cu":    f"{MOD}/_rec_het_cu.pdbqt",
     "hetatm_Zn":    f"{MOD}/_rec_het_zn.pdbqt",
     "no_metal":     f"{MOD}/hTYR_receptor_nocu.pdbqt"}
variant(V["hetatm_Cu"], "HETATM", "Cu")
variant(V["hetatm_Zn"], "HETATM", "Zn")

# probe: the PREPARED ligand (valid torsion tree), translated onto the metal
src = f"{LIG}/CHEMBL537.pdbqt"
lines = open(src).read().splitlines(True)
xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                for l in lines if l.startswith(("ATOM", "HETATM"))])
shift = mid - xyz.mean(0)
probe = f"{MOD}/_probe_on_metal.pdbqt"
with open(probe, "w") as f:
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            v = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) + shift
            f.write(l[:30] + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}" + l[54:])
        else:
            f.write(l)
ods = [min(np.linalg.norm(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                          + shift - c) for c in cu)
       for l in lines if l.startswith(("ATOM", "HETATM")) and l.split()[-1] in ("OA", "O")]
print(f"probe on metal: min O-Cu = {min(ods):.2f} A   ({len(ods)} oxygens)\n")

BASE = """-0.035579 gauss(o=0,_w=0.5,_c=8)
-0.005156 gauss(o=3,_w=2,_c=8)
0.840245 repulsion(o=0,_c=8)
-0.035069 hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439 non_dir_h_bond(g=-0.7,_b=0,_c=8)
"""
open(f"{MOD}/_plain.score", "w").write(BASE)
for t in ("Cu", "Zn"):
    open(f"{MOD}/_t_{t}.score", "w").write(
        BASE + f"-5.0 atom_type_gaussian(t1={t},t2=Oxygen,o=0,_w=1.5,_c=8)\n")

def score(rec, sf, show=False):
    p = subprocess.run(["smina", "-r", rec, "-l", probe, "--custom_scoring", sf,
                        "--score_only"], capture_output=True, text=True)
    if show:
        print("--- raw smina output ---")
        print((p.stdout or "")[-600:])
        if p.stderr.strip(): print("--- stderr ---\n" + p.stderr[-400:])
        print("------------------------\n")
    m = re.search(r"Affinity:\s*(-?\d+\.?\d*)", p.stdout)
    if m: return float(m.group(1)), ""
    err = (p.stderr or p.stdout).strip().splitlines()
    return None, (err[-1][:90] if err else "no output")

print("sanity check on one call:")
score(V["obabel_as_is"], f"{MOD}/_plain.score", show=True)

print(f"{'receptor':14s} {'base':>9s} {'+Cu-O':>9s} {'+Zn-O':>9s}   verdict")
for name, path in V.items():
    vals, errs = [], []
    for sf in ("_plain", "_t_Cu", "_t_Zn"):
        v, e = score(path, f"{MOD}/{sf}.score")
        vals.append(v); errs.append(e)
    if vals[0] is None:
        print(f"{name:14s} FAILED: {errs[0]}")
        continue
    seen = any(v is not None and abs(v - vals[0]) > 0.01 for v in vals[1:])
    fmt = lambda v: f"{v:9.3f}" if v is not None else f"{'err':>9s}"
    print(f"{name:14s} {fmt(vals[0])} {fmt(vals[1])} {fmt(vals[2])}   "
          f"{'METAL SEEN' if seen else 'metal invisible'}")
