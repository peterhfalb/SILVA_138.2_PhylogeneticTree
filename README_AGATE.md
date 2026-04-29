# SILVA 138.2 Tree Building - Agate Quick Start

Choose your path:

## Local Machine (Laptop/Desktop)

For testing or smaller datasets:

1. Follow [SETUP.md](SETUP.md)
2. Install Python packages and FastTree locally
3. Run: `python build_silva_tree.py SILVA_138.2_...fasta.gz`

## Agate Supercomputer (UMN MSI)

For the full SILVA 138.2 dataset (recommended):

1. Follow [AGATE_SETUP.md](AGATE_SETUP.md)
2. Upload this folder to Agate via SSH
3. Run setup: `bash slurm_setup.sh`
4. Submit job: `sbatch submit_tree_job.sbatch`

---

## File guide

| File | Purpose |
|------|---------|
| `build_silva_tree.py` | Main Python script (same for all environments) |
| `requirements.txt` | Python dependencies |
| `SETUP.md` | Local machine setup instructions |
| `slurm_setup.sh` | Agate one-time setup (installs venv + dependencies) |
| `submit_tree_job.sbatch` | SLURM job submission script |
| `AGATE_SETUP.md` | Complete Agate workflow guide |

---

## Typical Agate workflow

```bash
# 1. SSH to Agate
ssh agate.msi.umn.edu

# 2. Upload if needed or navigate to folder
cd ~/silva_tree_build

# 3. One-time setup (if first time)
bash slurm_setup.sh

# 4. Download SILVA alignment if needed
# NOTE: Must be the ALIGNED version (full_align), not unaligned sequences
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz

# 5. Edit submit_tree_job.sbatch (paths, resources)
nano submit_tree_job.sbatch

# 6. Submit to queue
sbatch submit_tree_job.sbatch

# 7. Monitor
squeue -u $USER
tail -f logs/silva-tree-*.out

# 8. Download results when done
# (from local machine)
scp -r agate.msi.umn.edu:~/silva_tree_build/results ~/
```

---

## Estimated times

| Dataset | Machine | Cores | Time |
|---------|---------|-------|------|
| Full SILVA 138.2 (510K seqs) | Agate | 16 | 1-2h |
| Full SILVA 138.2 | Laptop | 8 | 8-10h |
| Test subset (10K seqs) | Laptop | 4 | 5-10m |

---

**Next step**: Go to [AGATE_SETUP.md](AGATE_SETUP.md) to start!
