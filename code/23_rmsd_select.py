#!/usr/bin/env python3
"""Pick the metal term by symmetry-aware RMSD to the crystal pose, then
run attempt 4 on dev if (and only if) it beats baseline."""
import os, re, sys, json, math, subprocess, itertools
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD, LIG, RES, CUR = f"{ROOT}/models", f"{ROOT}/ligands", f"{ROOT}/results", f"{ROOT}/data/curated"
BASE = """-0.035579 gauss(o=0,_w=0.5,_c=8)
-0.005156 gauss(o=3,_w=2,_c=8)
0.840245 repulsion(o=0,_c=8)
-0.035069 hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439 non_dir_h_bond(g=-0.7,_b=0,_c=8)
"""
COORD = ["OxygenXSDonorAcceptor", "OxygenXSAcceptor",
         "NitrogenXSAcceptor", "NitrogenXSDonorAcceptor", "SulfurAcceptor"]
T2E = {"C": "C", "A": "C", "OA": "O", "O": "O", "N": "N", "NA": "N",
       "SA": "S", "S": "S", "F": "F", "Cl": "Cl", "Br": "Br", "I": "I"}

het = {}
for l in open(f"{MOD}/2Y9X.pdb"):
    if l.startswith("HETATM") and l[21] == "A":
        het.setdefault(l[17:20].strip(), []).append(l)
cry_lines = [l for l in het["0TR"] if l[76:78].strip() != "H"]
cry = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in cry_lines])
cry_el = [l[76:78].strip() for l in cry_lines]
ref_c = cry.mean(0)
cu = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in het["CU"]])[:2]

def sym_rmsd(pos, el):
    if len(pos) != len(cry): return float("nan")
    groups = {}
    for i, e in enumerate(el): groups.setdefault(e, []).append(i)
    cgroups = {}
    for i, e in enumerate(cry_el): cgroups.setdefault(e, []).append(i)
    if set(groups) != set(cgroups) or any(len(groups[k]) != len(cgroups[k]) for k in groups):
        return float("nan")
    keys = sorted(groups)
    if np.prod([math.factorial(len(groups[k])) for k in keys]) > 200000:
        return float(np.sqrt(((pos - cry) ** 2).sum(1).mean()))
    best = np.inf
    for combo in itertools.product(*[itertools.permutations(groups[k]) for k in keys]):
        m = np.empty(len(pos), dtype=int)
        for k, perm in zip(keys, combo):
            for ci, pi in zip(cgroups[k], perm): m[ci] = pi
        best = min(best, float(np.sqrt(((pos[m] - cry) ** 2).sum(1).mean())))
    return best

def redock(sf):
    out = f"{MOD}/_sel.pdbqt"
    if os.path.exists(out): os.remove(out)
    _p = subprocess.run(["smina", "-r", f"{MOD}/2Y9X_rec.pdbqt", "-l", f"{MOD}/0TR.pdbqt",
        "--center_x", f"{ref_c[0]:.2f}", "--center_y", f"{ref_c[1]:.2f}",
        "--center_z", f"{ref_c[2]:.2f}", "--size_x", "12", "--size_y", "12",
        "--size_z", "12", "--exhaustiveness", "16", "--num_modes", "1",
        "--seed", "42", "--cpu", "4", "-o", out]
        + (["--custom_scoring", sf] if sf else []), capture_output=True, text=True)
    if _p.returncode != 0 or not os.path.exists(out):
        raise SystemExit(f"smina FAILED (rc={_p.returncode}) scoring={sf}\n"
                         + (_p.stderr or _p.stdout)[-500:])
    pos, el, on = [], [], False
    for l in open(out):
        if l.startswith("MODEL"): on = True
        elif on and l.startswith(("ATOM", "HETATM")):
            t = l.split()[-1]
            if t == "HD": continue
            pos.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            el.append(T2E.get(t, t))
        elif l.startswith("ENDMDL") and pos: break
    return sym_rmsd(np.array(pos), el)

print(f"{'metal term':32s} {'RMSD to crystal':>16s}")
base_rmsd = redock(None)
print(f"{'(none - baseline)':32s} {base_rmsd:14.2f} A")
best = (base_rmsd, None)
for off, w, wt in itertools.product([-1.0, -0.5, 0.0], [0.6, 1.2], [-1.0, -3.0]):
    sf = f"{MOD}/_sel.score"
    open(sf, "w").write(BASE + "".join(
        f"{wt} atom_type_gaussian(t1=Cu,t2={t},o={off},_w={w},_c=8)\n" for t in COORD))
    if (off, w, wt) == (-1.0, 0.6, -1.0):
        print("--- first scoring file ---")
        print(open(sf).read().rstrip())
        print("--------------------------")
    r = redock(sf)
    print(f"o={off:+.1f} w={w:.1f} wt={wt:+.1f}{'':>12s} {r:14.2f} A")
    if r < best[0]: best = (r, (off, w, wt))

print(f"\nbaseline RMSD {base_rmsd:.2f} A   best {best[0]:.2f} A")
if best[1] is None or best[0] >= base_rmsd - 0.1:
    print("\nThe metal term does NOT improve pose recovery. Stopping here —")
    print("this is the stopping rule firing: move the ranking layer to a")
    print("ligand-based model on the curated dataset.")
    sys.exit(0)
if best[0] > 2.0:
    print(f"\nNote: best RMSD {best[0]:.2f} A is above the conventional 2.0 A")
    print("success threshold. Improvement is real but pose recovery is not solid.")

off, w, wt = best[1]
FROZEN = f"{MOD}/htyr_frozen.score"
open(FROZEN, "w").write(BASE + "".join(
    f"{wt} atom_type_gaussian(t1=Cu,t2={t},o={off},_w={w},_c=8)\n" for t in COORD))
print(f"\nFROZEN: o={off} w={w} weight={wt} -> {FROZEN}")
print("Selected on crystal RMSD alone. Potency was not consulted.\n")

box = json.load(open(f"{MOD}/box.json"))
dev = pd.read_csv(f"{CUR}/dev.csv")
rows = []
print(f"attempt 4: docking dev ({len(dev)}) with the frozen function")
for i, mid in enumerate(dev.molecule_chembl_id, 1):
    tag = mid.replace(":", "_"); lp = f"{LIG}/{tag}.pdbqt"
    if not os.path.exists(lp): continue
    op = f"{RES}/att4_{tag}_out.pdbqt"
    if not os.path.exists(op):
        subprocess.run(["smina", "-r", f"{MOD}/hTYR_receptor.pdbqt", "-l", lp,
            "--custom_scoring", FROZEN,
            "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
            "--center_z", str(box["center_z"]), "--size_x", "18", "--size_y", "18",
            "--size_z", "18", "--exhaustiveness", "8", "--num_modes", "5",
            "--seed", "42", "--cpu", "4", "-o", op], capture_output=True)
    sc = None
    for l in open(op):
        m = re.search(r"(?:minimizedAffinity|VINA RESULT:)\s+(-?\d+\.\d+)", l)
        if m: sc = float(m.group(1)); break
    if sc is not None: rows.append({"molecule_chembl_id": mid, "score": sc})
    if i % 10 == 0: print(f"  {i}/{len(dev)}", flush=True)

d = dev.merge(pd.DataFrame(rows), on="molecule_chembl_id")
d["heavy"] = d.smiles.map(lambda s: Chem.MolFromSmiles(s).GetNumHeavyAtoms()
                          if isinstance(s, str) and Chem.MolFromSmiles(s) else np.nan)
d.to_csv(f"{RES}/attempt4_dev.csv", index=False)
v = d.dropna(subset=["heavy", "pAct_median", "score"])
raw, p_raw = spearmanr(-v.score, v.pAct_median)
sl = linregress(v.heavy, v.score)
resid = v.score - (sl.slope * v.heavy + sl.intercept)
corr, p_corr = spearmanr(-resid, v.pAct_median)
k = max(1, int(round(0.20 * len(v))))
ef = (v.nsmallest(k, "score").pAct_median >= 5.5).mean() / max((v.pAct_median >= 5.5).mean(), 1e-9)
print(f"\n== attempt 4, dev (n={len(v)}) ==")
print(f"  raw rho            {raw:+.3f} (p={p_raw:.1e})")
print(f"  size-corrected rho {corr:+.3f} (p={p_corr:.1e})   <- gate condition 1'  (>= 0.30)")
print(f"  EF@20%             {ef:.2f}                       <- gate condition 2   (>= 2.0)")
