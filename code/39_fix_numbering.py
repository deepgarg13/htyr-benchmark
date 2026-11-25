#!/usr/bin/env python3
"""Recover the UniProt numbering by aligning the receptor's own sequence,
then redo Workstream 3 with residue features that mean what they say."""
import os, requests
import numpy as np, pandas as pd
from scipy.stats import spearmanr, linregress
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/hTYR_Platform")
MOD, RES, CUR = f"{ROOT}/models", f"{ROOT}/results", f"{ROOT}/data/curated"
THREE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
         "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
         "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

seq = requests.get("https://rest.uniprot.org/uniprotkb/P14679.fasta",
                   timeout=60).text.split("\n", 1)[1].replace("\n", "")
print(f"UniProt P14679: {len(seq)} residues")

obs = []
for l in open(f"{MOD}/hTYR_receptor.pdb"):
    if l.startswith("ATOM") and l[12:16].strip() == "CA":
        r = l[17:20].strip()
        if r in THREE: obs.append((int(l[22:26]), THREE[r]))
print(f"receptor: {len(obs)} residues, file numbering {obs[0][0]}-{obs[-1][0]}")

best = max(((sum(1 for n, a in obs
                 if 0 < n + off <= len(seq) and seq[n + off - 1] == a), off)
            for off in range(-400, 401)), key=lambda t: t[0])
hits, OFF = best
print(f"\nbest offset: uniprot = file {OFF:+d}   ({hits}/{len(obs)} residues match)")
if hits < 0.9 * len(obs):
    raise SystemExit("alignment failed - inspect the receptor PDB before continuing")

KEY = {180:"HIS",202:"HIS",211:"HIS",363:"HIS",367:"HIS",390:"HIS",
       377:"VAL",380:"SER",375:"SER",364:"ASN",368:"ILE",352:"THR",
       304:"HIS",306:"LYS",308:"ARG",347:"PHE",203:"GLU"}
fileno = {u: u - OFF for u in KEY}
found = {}
for l in open(f"{MOD}/hTYR_receptor.pdb"):
    if l.startswith("ATOM") and l[12:16].strip() == "CA":
        found[int(l[22:26])] = l[17:20].strip()
print("\n== key residues after correction ==")
ok = 0
for u in sorted(KEY):
    got = found.get(fileno[u], "ABSENT")
    good = got == KEY[u]; ok += good
    print(f"  {KEY[u]}{u:<4} -> file {fileno[u]:<4} {got:6s} {'ok' if good else 'MISMATCH'}")
print(f"  {ok}/{len(KEY)} correct")
if ok < len(KEY) - 1:
    raise SystemExit("key residues still wrong - stopping")

prot, cu = [], []
for l in open(f"{MOD}/hTYR_receptor.pdb"):
    if not l.startswith(("ATOM", "HETATM")): continue
    name, resn, rs = l[12:16].strip(), l[17:20].strip(), l[22:26].strip()
    xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    if resn == "CU": cu.append(xyz); continue
    if name.startswith("H") or l[76:78].strip() == "H": continue
    prot.append((int(rs) + OFF if rs.lstrip("-").isdigit() else -1, name[0], xyz))
cu = np.array(cu)[:2]; mid = cu.mean(0)
P = np.array([a[2] for a in prot]); P_res = np.array([a[0] for a in prot])
P_el = np.array([a[1] for a in prot])
KEYS = sorted(KEY)

def pose(path):
    xyz, el, on = [], [], False
    for l in open(path):
        if l.startswith("MODEL"): on = True
        elif l.startswith("ENDMDL") and xyz: break
        elif on and l.startswith(("ATOM", "HETATM")):
            t = l.split()[-1]
            if t == "HD": continue
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            el.append({"A":"C","OA":"O","NA":"N","SA":"S"}.get(t, t[0]))
    return (np.array(xyz), np.array(el)) if xyz else (None, None)

def feats(xyz, el):
    f = {}
    for i, tag in enumerate("AB"):
        d = np.linalg.norm(xyz - cu[i], axis=1)
        f[f"dCu{tag}_min"] = d.min()
        for e in "ONS":
            m = el == e
            f[f"dCu{tag}_{e}"] = d[m].min() if m.any() else 12.0
    f["depth"] = np.linalg.norm(xyz.mean(0) - mid)
    D = np.linalg.norm(xyz[:, None, :] - P[None, :, :], axis=2)
    close = D < 4.5
    pl = np.isin(el, ["O","N"])[:, None]; pp = np.isin(P_el, ["O","N"])[None, :]
    f["polar_contacts"] = int(((D < 3.5) & pl & pp).sum())
    n = len(xyz)
    for r in KEYS:                      # normalised: contacts per ligand atom
        f[f"res_{r}"] = close[:, P_res == r].sum() / n
    return f

df = pd.read_csv(f"{CUR}/multitask_train.csv")
hum = df[df.task == "human_enzyme"].reset_index(drop=True)
rows, keep = [], []
for i, cid in enumerate(hum.cid):
    p = f"{RES}/att4_{cid.replace(':','_')}_out.pdbqt"
    if not os.path.exists(p): continue
    xyz, el = pose(p)
    if xyz is None: continue
    rows.append(feats(xyz, el)); keep.append(i)
hum = hum.iloc[keep].reset_index(drop=True)
IFP = pd.DataFrame(rows)
print(f"\nfingerprints: {len(IFP)} compounds x {IFP.shape[1]} features")

DESC = ["MolWt","MolLogP","TPSA","NumHDonors","NumHAcceptors","NumRotatableBonds",
        "NumAromaticRings","FractionCSP3","HeavyAtomCount","RingCount","NHOHCount","NOCount"]
fn = {n: f for n, f in Descriptors.descList if n in DESC}
Xd = np.array([[fn[n](Chem.MolFromSmiles(s)) for n in DESC] for s in hum.smiles])
y = hum.pAct.values; heavy = Xd[:, DESC.index("HeavyAtomCount")]
Xr = IFP[[c for c in IFP.columns if c.startswith("res_")]].values
Xm = IFP[[c for c in IFP.columns if c.startswith("dCu")] + ["depth","polar_contacts"]].values

def score(X, seed):
    oof = np.zeros(len(y))
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        m = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1).fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    sl = linregress(heavy, oof)
    return spearmanr(oof - (sl.slope*heavy + sl.intercept), y).statistic

print(f"\n== dev CV, n={len(y)}, 20 seeds ==")
for label, X in (("A: descriptors only", Xd),
                 ("B: descriptors + residue contacts", np.hstack([Xd, Xr])),
                 ("C: residue contacts alone", Xr),
                 ("D: descriptors + all IFP", np.hstack([Xd, Xr, Xm]))):
    s = [score(X, 42 + r) for r in range(20)]
    print(f"  {label:36s} rho = {np.mean(s):+.3f} +/- {np.std(s):.3f}")

m = RandomForestRegressor(n_estimators=600, min_samples_leaf=2,
                          random_state=42, n_jobs=-1).fit(Xr, y)
print("\n== which residues carry signal? ==")
SPEC = {377,380,375,364,368,352,304,306,308}
for n, imp in sorted(zip([c for c in IFP.columns if c.startswith("res_")],
                         m.feature_importances_), key=lambda t: -t[1])[:8]:
    u = int(n[4:])
    print(f"  {KEY[u]}{u:<5} {imp:.3f}{'   <- hsTYR-specific' if u in SPEC else ''}")
