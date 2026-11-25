#!/usr/bin/env bash
echo "== docking engines =="
for b in vina smina gnina qvina qvina2 autodock4 autodock_gpu autogrid4 rbdock; do
  p=$(command -v $b 2>/dev/null) && echo "  $b -> $p" || echo "  $b -> MISSING"
done

echo -e "\n== chemistry / structure CLI =="
for b in obabel pdb2pqr reduce pymol mk_prepare_ligand.py prepare_receptor4.py; do
  p=$(command -v $b 2>/dev/null) && echo "  $b -> $p" || echo "  $b -> MISSING"
done

echo -e "\n== AmberTools =="
for b in antechamber parmchk2 tleap sander pmemd; do
  p=$(command -v $b 2>/dev/null) && echo "  $b -> $p" || echo "  $b -> MISSING"
done

echo -e "\n== python packages (current env: $CONDA_DEFAULT_ENV) =="
python - <<'PY'
import importlib.metadata as md
for p in ["rdkit","openmm","openff-toolkit","pdbfixer","biopython","mdanalysis",
          "mdtraj","prolif","meeko","vina","scipy","scikit-learn","torch",
          "torch-geometric","py3Dmol","gemmi","parmed"]:
    try: print(f"  {p:20s} {md.version(p)}")
    except md.PackageNotFoundError: print(f"  {p:20s} MISSING")
PY

echo -e "\n== envs =="
conda env list
echo -e "\n== gpu =="
command -v nvidia-smi >/dev/null && nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || echo "  no nvidia-smi (CPU only)"
