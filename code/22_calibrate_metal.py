#!/usr/bin/env python3
"""Attempt 4, step 1: find N/S type names, then calibrate a Cu-coordination
term against the 2Y9X CRYSTAL POSE — never against potency."""
import os, re, subprocess, itertools
import numpy as np

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD = f"{ROOT}/models"
BASE = """-0.035579 gauss(o=0,_w=0.5,_c=8)
-0.005156 gauss(o=3,_w=2,_c=8)
0.840245 repulsion(o=0,_c=8)
-0.035069 hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439 non_dir_h_bond(g=-0.7,_b=0,_c=8)
"""

print("== type names for N and S acceptors ==")
CAND = ["Nitrogen", "NitrogenXSAcceptor", "NitrogenXSDonorAcceptor",
        "NitrogenXSDonor", "Sulfur", "SulfurAcceptor", "NA", "SA", "S"]
def last_term(sf, rec, lig):
    p = subprocess.run(["smina", "-r", rec, "-l", lig, "--custom_scoring", sf,
                        "--score_only"], capture_output=True, text=True)
    m = re.search(r"Term values, before weighting:\s*\n##\s*(.+)", p.stdout)
    return float(m.group(1).split()[-1]) if m else None

good = ["OxygenXSDonorAcceptor"]
for t2 in CAND:
    sf = f"{MOD}/_p.score"
    open(sf, "w").write(BASE + f"-1.0 atom_type_gaussian(t1=Cu,t2={t2},o=0,_w=1.5,_c=8)\n")
    v = last_term(sf, f"{MOD}/hTYR_receptor.pdbqt", f"{MOD}/_probe_on_metal.pdbqt")
    status = "accepted" if v is not None else "REJECTED"
    print(f"  t2={t2:26s} {status}")
    if v is not None and t2 not in ("NA", "SA", "S"):
        good.append(t2)
coord_types = ["OxygenXSDonorAcceptor", "OxygenXSAcceptor",
               "NitrogenXSAcceptor", "NitrogenXSDonorAcceptor", "SulfurAcceptor"]
coord_types = [t for t in coord_types if t in good or t == "OxygenXSDonorAcceptor"]
print(f"\n  coordinating types to include: {coord_types}")

# --- crystal reference ---
het = {}
for l in open(f"{MOD}/2Y9X.pdb"):
    if l.startswith("HETATM") and l[21] == "A":
        het.setdefault(l[17:20].strip(), []).append(l)
xyz_of = lambda ls: np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in ls])
cry = xyz_of(het["0TR"]); cu = xyz_of(het["CU"])[:2]
cry_o = xyz_of([l for l in het["0TR"] if l[76:78].strip() == "O"])
ref_c = cry.mean(0)
ref_ocu = min(np.linalg.norm(o - m) for o in cry_o for m in cu)
print(f"\n== crystal reference: centroid {ref_c.round(2)}  min O-Cu {ref_ocu:.2f} A ==")

def redock(sf):
    out = f"{MOD}/_cal.pdbqt"
    subprocess.run(["smina", "-r", f"{MOD}/2Y9X_rec.pdbqt", "-l", f"{MOD}/0TR.pdbqt",
        "--center_x", f"{ref_c[0]:.2f}", "--center_y", f"{ref_c[1]:.2f}",
        "--center_z", f"{ref_c[2]:.2f}", "--size_x", "12", "--size_y", "12",
        "--size_z", "12", "--exhaustiveness", "16", "--num_modes", "1",
        "--seed", "42", "--cpu", "4", "-o", out]
        + (["--custom_scoring", sf] if sf else []), capture_output=True, text=True)
    pos, o, on = [], [], False
    for l in open(out):
        if l.startswith("MODEL"): on = True
        elif on and l.startswith(("ATOM", "HETATM")):
            v = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            pos.append(v)
            if l.split()[-1] in ("OA", "O"): o.append(v)
        elif l.startswith("ENDMDL") and pos: break
    pos = np.array(pos)
    ocu = min(np.linalg.norm(np.array(a) - b) for a in o for b in cu) if o else float("nan")
    return np.linalg.norm(pos.mean(0) - ref_c), ocu

print(f"\n{'metal term':34s} {'centroid off':>13s} {'min O-Cu':>10s}  {'err vs crystal':>14s}")
d, ocu = redock(None)
print(f"{'(none - baseline)':34s} {d:12.2f} A {ocu:9.2f} A {abs(ocu-ref_ocu):13.2f}")
results = []
for off, w, wt in itertools.product([-1.0, -0.5, 0.0], [0.6, 1.2], [-1.0, -3.0]):
    sf = f"{MOD}/_cal.score"
    terms = "".join(f"{wt} atom_type_gaussian(t1=Cu,t2={t},o={off},_w={w},_c=8)\n"
                    for t in coord_types)
    open(sf, "w").write(BASE + terms)
    d, ocu = redock(sf)
    err = abs(ocu - ref_ocu) + abs(d)
    results.append((err, off, w, wt, d, ocu))
    print(f"o={off:+.1f} w={w:.1f} wt={wt:+.1f}{'':>14s} {d:12.2f} A {ocu:9.2f} A {err:13.2f}")

best = min(results)
print(f"\nBEST by structural agreement: o={best[1]} w={best[2]} weight={best[3]}")
print(f"  centroid off {best[4]:.2f} A, min O-Cu {best[5]:.2f} A vs crystal {ref_ocu:.2f} A")
sf = f"{MOD}/htyr_metal.score"
open(sf, "w").write(BASE + "".join(
    f"{best[3]} atom_type_gaussian(t1=Cu,t2={t},o={best[1]},_w={best[2]},_c=8)\n"
    for t in coord_types))
print(f"\nfrozen scoring function -> {sf}")
print("Calibrated on crystal geometry only. Potency data was not consulted.")
