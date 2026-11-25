#!/usr/bin/env python3
"""Transplant the binuclear Cu site from 2Y9X onto the AlphaFold model of
human TYR, by superposing the conserved six-histidine NE2 cage."""
import os, itertools, requests
import numpy as np

MOD = os.path.expanduser("~/hTYR_Platform/models")
AF, TPL = f"{MOD}/AF_P14679.pdb", f"{MOD}/2Y9X.pdb"

TPL_CUA = [("A", 61), ("A", 85), ("A", 94)]
TPL_CUB = [("A", 259), ("A", 263), ("A", 296)]

print("== UniProt copper-binding residues for P14679 ==")
j = requests.get("https://rest.uniprot.org/uniprotkb/P14679.json", timeout=60).json()
import re
BIND = [f for f in j.get("features", []) if f.get("type") == "Binding site"]
print(f"  {len(BIND)} binding-site features in total:")
for f in BIND:
    lg = f.get("ligand") or {}
    print(f"    res {f['location']['start']['value']:4d}  name={lg.get('name')!r}"
          f"  id={lg.get('id')!r}  label={lg.get('label')!r}")

CU = re.compile(r"\bcu\b|copper|chebi:29036|chebi:49552", re.I)
his = []
for f in BIND:
    lg = f.get("ligand") or {}
    blob = " ".join(str(v) for v in (lg.get("name"), lg.get("id"), lg.get("label"),
                                     f.get("description")) if v)
    if CU.search(blob):
        his.append((f["location"]["start"]["value"], lg.get("name"),
                    f.get("description", "")))
print(f"\n  matched {len(his)} copper-binding residues:")
for pos, lig, desc in his:
    print(f"    residue {pos:4d}  {lig}  {desc}")
nums = sorted({p for p, _, _ in his})
if len(nums) != 6:
    raise SystemExit(f"\nSTOP: need exactly 6 copper-binding residues, got "
                     f"{len(nums)}: {nums}\nCheck the feature dump above.")
af_cua, af_cub = nums[:3], nums[3:]
print(f"\n  CuA triad (human): {af_cua}\n  CuB triad (human): {af_cub}")

def atoms(path):
    out = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            out.append((l[:6].strip(), l[12:16].strip(), l[17:20].strip(), l[21],
                        int(l[22:26]), np.array([float(l[30:38]), float(l[38:46]),
                                                 float(l[46:54])]), l))
    return out

def ne2(ats, chain, resseq):
    for rec, name, resn, ch, rs, xyz, _ in ats:
        if rec == "ATOM" and name == "NE2" and resn == "HIS" and ch == chain and rs == resseq:
            return xyz
    raise SystemExit(f"NE2 not found for HIS {chain}{resseq}")

def plddt(ats, resseq):
    for rec, name, resn, ch, rs, xyz, line in ats:
        if rec == "ATOM" and name == "CA" and rs == resseq:
            return float(line[60:66])
    return float("nan")

a_af, a_tp = atoms(AF), atoms(TPL)

print("\n== AlphaFold confidence at the active site ==")
for n in nums:
    print(f"  residue {n:4d}  pLDDT {plddt(a_af, n):5.1f}")

cu_tpl = [xyz for rec, name, resn, ch, rs, xyz, _ in a_tp
          if rec == "HETATM" and resn == "CU" and ch == "A"]
print(f"\ntemplate Cu pair, Cu–Cu = {np.linalg.norm(cu_tpl[0]-cu_tpl[1]):.2f} A")

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Qc - R @ Pc

tplA = [ne2(a_tp, c, r) for c, r in TPL_CUA]
tplB = [ne2(a_tp, c, r) for c, r in TPL_CUB]
afA  = [ne2(a_af, "A", r) for r in af_cua]
afB  = [ne2(a_af, "A", r) for r in af_cub]

best = None
for pa in itertools.permutations(range(3)):
    for pb in itertools.permutations(range(3)):
        P = np.array([tplA[i] for i in pa] + [tplB[i] for i in pb])
        Q = np.array(afA + afB)
        R, t = kabsch(P, Q)
        rmsd = np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean())
        if best is None or rmsd < best[0]:
            best = (rmsd, R, t, pa, pb)

rmsd, R, t, pa, pb = best
print(f"\nbest six-point superposition RMSD: {rmsd:.2f} A   pairing {pa}/{pb}")
if rmsd > 1.5:
    print("  !! high RMSD — the His cage may not correspond; inspect before docking")

cu_new = [R @ c + t for c in cu_tpl]
print(f"placed Cu–Cu distance: {np.linalg.norm(cu_new[0]-cu_new[1]):.2f} A")
print("\n== validation: Cu to NE2 distances (expect ~2.0–2.3 A) ==")
for label, cu, triad in (("CuA", cu_new[0], afA), ("CuB", cu_new[1], afB)):
    ds = [np.linalg.norm(cu - x) for x in triad]
    flag = "" if all(1.7 < d < 2.8 for d in ds) else "   <-- CHECK"
    print(f"  {label}: " + ", ".join(f"{d:.2f}" for d in ds) + flag)

out = f"{MOD}/hTYR_Cu.pdb"
with open(out, "w") as f:
    ser = 0
    for l in open(AF):
        if l.startswith("ATOM"):
            ser = int(l[6:11]); f.write(l)
    for i, cu in enumerate(cu_new):
        ser += 1
        f.write(f"HETATM{ser:5d} CU    CU A{900+i:4d}    "
                f"{cu[0]:8.3f}{cu[1]:8.3f}{cu[2]:8.3f}  1.00  0.00          CU\n")
    f.write("END\n")
print(f"\nwrote {out}")
