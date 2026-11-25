#!/usr/bin/env python3
"""Dock the locked benchmark under three receptor/scoring conditions."""
import os, re, glob, subprocess, json, sys
import pandas as pd

MOD = os.path.expanduser("~/hTYR_Platform/models")
LIG = os.path.expanduser("~/hTYR_Platform/ligands")
RES = os.path.expanduser("~/hTYR_Platform/results")
CUR = os.path.expanduser("~/hTYR_Platform/data/curated")
os.makedirs(RES, exist_ok=True)
box = json.load(open(f"{MOD}/box.json"))

# receptor variant with formal Cu charge
src, chg = f"{MOD}/hTYR_receptor.pdbqt", f"{MOD}/hTYR_receptor_cu2.pdbqt"
if not os.path.exists(chg):
    with open(chg, "w") as f:
        for l in open(src):
            if l.startswith(("ATOM", "HETATM")) and l.rstrip().endswith("Cu"):
                l = l[:70] + f"{2.000:7.3f}" + l[77:]
            f.write(l)
    print(f"wrote {chg} (Cu q=+2.00)")

RUNS = [("A_vina",   src, "vina"),
        ("B_ad4_q0", src, "ad4_scoring"),
        ("C_ad4_q2", chg, "ad4_scoring")]

EXH = int(sys.argv[1]) if len(sys.argv) > 1 else 8
ligs = sorted(glob.glob(f"{LIG}/*.pdbqt"))
print(f"{len(ligs)} ligands, exhaustiveness {EXH}")

for name, rec, scoring in RUNS:
    out_csv = f"{RES}/scores_{name}.csv"
    done = set()
    if os.path.exists(out_csv):
        done = set(pd.read_csv(out_csv).molecule_chembl_id)
    rows = []
    print(f"\n== run {name} ({scoring}) ==", flush=True)
    for i, lp in enumerate(ligs, 1):
        tag = os.path.basename(lp)[:-6]
        mid = tag.replace("PUBCHEM_", "PUBCHEM:")
        if mid in done:
            continue
        cmd = ["smina", "-r", rec, "-l", lp, "--scoring", scoring,
               "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
               "--center_z", str(box["center_z"]), "--size_x", str(box["size_x"]),
               "--size_y", str(box["size_y"]), "--size_z", str(box["size_z"]),
               "--exhaustiveness", str(EXH), "--num_modes", "5",
               "--seed", "42", "--cpu", "4", "-o", f"{RES}/{name}_{tag}_out.pdbqt"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        m = re.search(r"^\s*1\s+(-?\d+\.\d+)", p.stdout, re.M)
        if m:
            rows.append({"molecule_chembl_id": mid, "score": float(m.group(1))})
        else:
            print(f"  {tag}: no score\n{p.stderr[:200]}")
        if i % 10 == 0:
            print(f"  {i}/{len(ligs)}", flush=True)
    if rows:
        new = pd.DataFrame(rows)
        if os.path.exists(out_csv):
            new = pd.concat([pd.read_csv(out_csv), new], ignore_index=True)
        new.to_csv(out_csv, index=False)
    print(f"  wrote {out_csv}")

sys.path.insert(0, os.path.expanduser("~/hTYR_Platform/code"))
from gate import run_gate
bench = pd.read_csv(f"{CUR}/benchmark_v2.csv")
for name, _, _ in RUNS:
    f = f"{RES}/scores_{name}.csv"
    if not os.path.exists(f):
        continue
    print(f"\n{'='*58}\nGATE — {name}\n{'='*58}")
    run_gate(bench, pd.read_csv(f), score_col="score", higher_is_better=False)
