#!/usr/bin/env python3
"""Find the atom-type names that actually match our Cu and ligand oxygens,
by reading smina's own per-term values."""
import os, re, subprocess, itertools

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD = f"{ROOT}/models"
REC = f"{MOD}/hTYR_receptor.pdbqt"
PROBE = f"{MOD}/_probe_on_metal.pdbqt"

BASE = """-0.035579 gauss(o=0,_w=0.5,_c=8)
-0.005156 gauss(o=3,_w=2,_c=8)
0.840245 repulsion(o=0,_c=8)
-0.035069 hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439 non_dir_h_bond(g=-0.7,_b=0,_c=8)
"""

METALS = ["Cu", "Zn", "Zinc", "GenericMetal", "Magnesium", "Manganese",
          "Calcium", "Iron", "Boron"]
OXY    = ["Oxygen", "OxygenXSAcceptor", "OxygenXSDonorAcceptor",
          "OxygenXSDonor", "OA", "O"]

def last_term(sf):
    p = subprocess.run(["smina", "-r", REC, "-l", PROBE, "--custom_scoring", sf,
                        "--score_only"], capture_output=True, text=True)
    m = re.search(r"Term values, before weighting:\s*\n##\s*(.+)", p.stdout)
    if not m:
        err = (p.stderr or p.stdout).strip().splitlines()
        return None, (err[-1][:70] if err else "no output")
    vals = m.group(1).split()
    return float(vals[-1]), ""

print("searching for type names that match (non-zero last term = MATCH)\n")
hits = []
for t1, t2 in itertools.product(METALS, OXY):
    sf = f"{MOD}/_probe.score"
    open(sf, "w").write(
        BASE + f"-1.0 atom_type_gaussian(t1={t1},t2={t2},o=0,_w=1.5,_c=8)\n")
    v, err = last_term(sf)
    if v is None:
        continue
    if abs(v) > 1e-6:
        print(f"  MATCH   t1={t1:14s} t2={t2:22s} raw term = {v:.5f}")
        hits.append((t1, t2, v))

if not hits:
    print("  no combination matched.")
    print("\n  falling back: what types does smina assign? listing receptor atom")
    print("  types near the metal via --atom_terms")
    p = subprocess.run(["smina", "-r", REC, "-l", PROBE, "--score_only",
                        "--atom_terms", f"{MOD}/_atom_terms.txt"],
                       capture_output=True, text=True)
    if os.path.exists(f"{MOD}/_atom_terms.txt"):
        print("".join(open(f"{MOD}/_atom_terms.txt").readlines()[:25]))
    else:
        print("  (--atom_terms produced nothing)")
        print((p.stderr or p.stdout)[-300:])
else:
    best = max(hits, key=lambda h: abs(h[2]))
    print(f"\n  use: t1={best[0]} t2={best[1]}")
