# Running SILVA Tree Building on Agate (UMN MSI)

Complete guide for building SILVA 138.2 phylogenetic trees on the University of Minnesota's Agate supercomputer.

## Overview

This setup allows you to submit a SLURM job that will:
1. Load required software modules
2. Install Python dependencies in a virtual environment
3. Build the phylogenetic tree with FastTree
4. Save results for download and visualization

## Prerequisites

- SSH access to Agate: `falb0011@agate.msi.umn.edu`
- A valid allocation on Agate (check with `salloc --help`)
- The SILVA 138.2 alignment file downloaded or accessible on Agate

## Step 1: Upload files to Agate

From your local machine:

```bash
# Copy entire project folder to Agate
scp -r /path/to/silva_tree_build agate.msi.umn.edu:~/

# Or if already uploaded, SSH in
ssh agate.msi.umn.edu
cd ~/silva_tree_build
```

## Step 2: One-time setup

Run the setup script via SSH (only needed once):

```bash
# SSH into Agate
ssh agate.msi.umn.edu

# Navigate to project directory
cd ~/silva_tree_build

# Make setup script executable and run it
chmod +x slurm_setup.sh
bash slurm_setup.sh
```

This will:
- Load Python module
- Create a virtual environment (`venv/`)
- Install scikit-bio and numpy
- Check for FastTree availability

## Step 3: Download SILVA 138.2 alignment

If not already on Agate:

```bash
# Still in ~/silva_tree_build
# IMPORTANT: Make sure to download the ALIGNED version (with "full_align" in the name)
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz

# Verify download (~1.4 GB compressed, ~25 GB uncompressed)
ls -lh SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz
gunzip -t SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz && echo "✓ Archive is valid"
```

⚠️ **Important**: Make sure you're downloading the `full_align` version, not the unaligned sequences. The filename must contain `full_align`.

Alternatively, use your own file and update the path in `submit_tree_job.sbatch`.

## Step 4: Configure and submit SLURM job

Before submitting, edit `submit_tree_job.sbatch`:

```bash
# Edit configuration section
nano submit_tree_job.sbatch
```

Key parameters to adjust:

```bash
# Line 20: Path to SILVA alignment (if different)
SILVA_FILE="/home/$USER/silva_tree_build/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz"

# Line 23: Output directory for results
OUTPUT_DIR="/home/$USER/silva_tree_build/results"

# Line 32: Job name (optional)
#SBATCH --job-name=silva-tree

# Line 36: Cores to use (more = faster, longer queue wait)
# Recommended: 16 cores = ~1-2 hours runtime
#SBATCH --cpus-per-task=16

# Line 40: Memory (need ~12-16 GB for full dataset)
#SBATCH --mem=20G

# Line 44: Time limit (3 hours is safe default)
#SBATCH --time=03:00:00
```

### Resource recommendations

| Cores | Time | Memory | Queue Wait |
|-------|------|--------|-----------|
| 4     | 4-5h | 16G    | Short     |
| 8     | 2-3h | 16G    | Short     |
| 16    | 1-2h | 20G    | Medium    |
| 32    | 45m  | 24G    | Long      |

Submit the job:

```bash
# Make script executable
chmod +x submit_tree_job.sbatch

# Submit to queue
sbatch submit_tree_job.sbatch

# You'll see: "Submitted batch job <JOB_ID>"
```

## Step 5: Monitor job progress

```bash
# Check job status
squeue -u $USER

# View job details
scontrol show job <JOB_ID>

# Stream output in real-time
tail -f logs/silva-tree-*.out

# Check for errors
tail -f logs/silva-tree-*.err

# When complete, job will disappear from squeue
```

## Step 6: Download and visualize results

Once complete, download results to your local machine:

```bash
# From local machine (not Agate)
scp -r agate.msi.umn.edu:~/silva_tree_build/results ~/silva_results
```

Results folder contains:

- `silva-138.2-rooted-tree.nwk` - **Main output** (Newick format)
- `silva-138.2-tree.nwk` - Unrooted version
- `silva-138.2-masked-aln.fasta` - Masked alignment (large, optional)
- Logs with performance metrics

### Visualize the tree

Upload `silva-138.2-rooted-tree.nwk` to:
- **iTOL** (web): https://itol.embl.de/
  - Upload tree file
  - Explore interactive visualization
  - Export high-quality figures for publications

Or use desktop software:
- **FigTree**: http://tree.bio.ed.ac.uk/software/figtree/
- **Dendroscope**: http://dendroscope.org/

## Troubleshooting

### Job fails with "FastTree: command not found"

SLURM script includes `module load fasttree`, but if not available:

```bash
# Check available modules
module avail fasttree

# Load if available with different name
module load fasttree/2.1.10  # or similar version

# If not available, email MSI support
```

### "Out of memory" error

Increase memory in `submit_tree_job.sbatch`:

```bash
#SBATCH --mem=32G  # Increase from 20G to 32G
```

Then resubmit.

### Job times out

Increase time limit in `submit_tree_job.sbatch`:

```bash
#SBATCH --time=04:00:00  # Increase from 03:00:00 to 4 hours
```

Also try reducing cores (more cores = more memory contention).

### "Cannot access file" errors

Verify file paths in `submit_tree_job.sbatch`:

```bash
# SSH to Agate and check
ls -lh ~/silva_tree_build/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
```

Path must be absolute (start with `/home/`) not relative.

## Advanced usage

### Test with smaller dataset first

Create a small test alignment:

```bash
# Extract first 10K sequences for quick test
zcat SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz | head -40000 > silva_test.fasta

# Test locally
python3 build_silva_tree.py silva_test.fasta -t 4

# Once verified, run on full dataset via SLURM
```

### Run multiple analyses in parallel

Create multiple SLURM scripts with different parameters:

```bash
# submit_tree_job_16cores.sbatch
sbatch submit_tree_job_16cores.sbatch

# submit_tree_job_32cores.sbatch
sbatch submit_tree_job_32cores.sbatch

# Monitor both
squeue -u $USER
```

### Customize gap threshold

Edit `submit_tree_job.sbatch` and add to `python3 build_silva_tree.py` call:

```bash
python3 build_silva_tree.py \
    "$SILVA_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --gap-threshold 99.0  # Change from 99.56 to 99.0 for stricter filtering
```

## Performance benchmarks

Observed runtimes on Agate (full SILVA 138.2, 510K sequences):

| Cores | Wall Time | Est. Cost |
|-------|-----------|-----------|
| 4     | 5.5 hours | ~0.25 SU* |
| 8     | 2.5 hours | ~0.15 SU  |
| 16    | 1.5 hours | ~0.18 SU  |
| 32    | 45 min    | ~0.22 SU  |

*Service Units - check your allocation with `atq` or MSI dashboard

## Getting help

For issues with:
- **Tree building script**: Check logs in `logs/` directory
- **Agate/SLURM**: Email MSI support: hpc@umn.edu
- **SILVA data**: https://www.arb-silva.de/
- **FastTree**: http://microbesonline.org/fasttree/

## References

- Agate documentation: https://www.msi.umn.edu/content/agate
- SLURM documentation: https://slurm.schedmd.com/
- SILVA database: https://www.arb-silva.de/
- FastTree: http://microbesonline.org/fasttree/
