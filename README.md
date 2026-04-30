# SILVA 138.2 NR99 Phylogenetic Tree Pipeline

Build a high-quality phylogenetic tree from the SILVA NR99 138.2 rRNA database using modern sequence analysis tools.

## Overview

This pipeline constructs a reference phylogenetic tree by:

1. **Downloading** the unaligned SILVA NR99 138.2 sequences (~510K sequences)
2. **Clustering** at 97% sequence identity using VSEARCH to reduce redundancy
3. **Aligning** centroid sequences to the RF00177 SSU rRNA covariance model using cmalign (Infernal)
4. **Masking** high-gap columns in the alignment (keeping columns with <99.56% gaps)
5. **Building** a maximum-likelihood tree with IQ-TREE, using:
   - ModelFinder to select the best evolutionary model
   - 1000 ultrafast bootstrap replicates
6. **Pruning** the tree to 10,000 tips using Treemmer, maximizing phylogenetic diversity

All steps run on Agate (UMN HPC) as modular SLURM jobs.

### Output Files

| File | Description |
|------|-------------|
| `results/centroids_masked.fasta` | Final masked alignment for tree building |
| `results/iqtree/silva_138.2_nr99.treefile` | Maximum-likelihood tree (Newick format) |
| `results/iqtree/silva_138.2_nr99.log` | IQ-TREE detailed log |
| `results/treemmer/best/silva_138.2_nr99_10k.nwk` | Pruned tree (10,000 tips) |
| `results/centroids_taxonomy.tsv` | Taxonomic assignment per centroid sequence |
| `logs/pipeline.log` | Summary log with statistics from each step |

## Setup (once, interactive)

```bash
# On Agate login node:
cd /path/to/SILVA_138.2_NR99_fastTree
bash scripts/01_download_silva.sh
```

This downloads the ~2.2 GB SILVA FASTA file. If it already exists, the script verifies it.

## Installation

```bash
# Load Python and create virtual environment
module load python3
python3 -m venv venv
source venv/bin/activate
pip install -U pip
pip install -r requirements.txt
deactivate
```

## Running the Pipeline

### Option 1: Submit all jobs with dependencies (recommended)

```bash
bash scripts/submit_pipeline.sh
```

This submits all 4 jobs with dependencies, so they run sequentially. Each job waits for the previous one to complete successfully.

**Monitor progress:**
```bash
squeue -u $USER
slog <job_id>  # Stream logs for a specific job
```

### Option 2: Submit each step manually

```bash
# Step 2: Clustering
sbatch scripts/02_cluster.sbatch

# Step 3: Alignment & Masking (only after step 2 completes)
sbatch scripts/03_align_and_mask.sbatch

# Step 4: Tree Building (only after step 3 completes)
sbatch scripts/04_build_tree.sbatch

# Step 5: Pruning (only after step 4 completes)
sbatch scripts/05_prune_tree.sbatch
```

## Pipeline Details

### Step 2: Clustering (VSEARCH)

**File:** `scripts/02_cluster.sbatch`  
**Resources:** agsmall, 32 cores, 128 GB, 4 hours  
**Input:** SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz (~510K sequences)  
**Output:** 
- `results/centroids_97.fasta` — representative sequences at 97% identity
- `results/clusters_97.uc` — cluster assignments

Reduces ~510K sequences to ~50K–150K centroids (typical ~20% compression).

### Step 3: Alignment & Masking

**File:** `scripts/03_align_and_mask.sbatch`  
**Resources:** aglarge, 64 cores, 256 GB, 12 hours

**Substeps:**
1. **cmalign** — aligns centroid sequences to RF00177.cm (SSU rRNA bacteria model)
   - `--outformat afa` outputs aligned FASTA
   - `--cpu 64` uses all available cores
   - Output: `results/centroids_aligned.fasta` (~1500 columns)

2. **Masking** — removes high-gap columns using `mask_alignment.py`
   - Keeps columns with gap% < 99.56%
   - Output: `results/centroids_masked.fasta` (~1000–1200 columns after masking)

3. **Taxonomy Extraction** — parses SILVA headers to extract full lineage
   - Output: `results/centroids_taxonomy.tsv` (seqid + full semicolon-delimited lineage)
   - Fixes the bug from the original pipeline (extracts all ranks, not just last)

### Step 4: Tree Building (IQ-TREE)

**File:** `scripts/04_build_tree.sbatch`  
**Resources:** aglarge, 64 cores, 500 GB, 72 hours

**Key parameters:**
- `-m TEST` — runs ModelFinder to select best substitution model
- `-bb 1000` — 1000 ultrafast bootstrap replicates (fast, still statistically sound)
- `-nt AUTO -ntmax 64` — auto-detect optimal thread count

**Output:**
- `results/iqtree/silva_138.2_nr99.treefile` — ML tree (primary output)
- `results/iqtree/silva_138.2_nr99.log` — run log
- `results/iqtree/silva_138.2_nr99.iqtree` — detailed report with selected model

Expected runtime: 48–72 hours depending on alignment size.

### Step 5: Tree Pruning (Treemmer)

**File:** `scripts/05_prune_tree.sbatch`  
**Resources:** agsmall, 8 cores, 32 GB, 4 hours

**Parallelism:** Runs 4 independent Treemmer instances in parallel (using GNU parallel), each with a different random seed. The best result is selected.

**Algorithm:** Iteratively removes tips that contribute least to phylogenetic diversity (PD), preserving maximum diversity in the final 10K-tip tree.

**Output:**
- `results/treemmer/best/silva_138.2_nr99_10k.nwk` — pruned tree (10,000 tips)
- `results/treemmer/best/silva_138.2_nr99_10k_tips.txt` — list of retained tip labels

## Troubleshooting

### Job fails in Step 2 (VSEARCH)
- Check: `logs/silva-cluster-*.err`
- Ensure SILVA FASTA is downloaded: `ls -lh SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz`

### Job fails in Step 3 (cmalign)
- Check: `logs/silva-align-*.err`
- cmalign memory usage can be high; aglarge partition ensures sufficient RAM
- If memory error: increase `--mem` in the sbatch script

### Job fails in Step 4 (IQ-TREE)
- Check: `logs/silva-tree-*.err` and `results/iqtree/silva_138.2_nr99.log`
- IQ-TREE is memory-hungry on large alignments; 500 GB is usually sufficient
- Run time: 48–72 hours for ~100K sequences is expected

### Job fails in Step 5 (Treemmer)
- Check: `logs/silva-prune-*.err`
- Treemmer Python package might need reinstalling: `pip install -U treemmer`
- Ensure tree file was created in step 4: `ls -l results/iqtree/silva_138.2_nr99.treefile`

## Repository Structure

```
scripts/
  ├── 01_download_silva.sh        # Download SILVA FASTA
  ├── 02_cluster.sbatch           # VSEARCH clustering
  ├── 03_align_and_mask.sbatch    # cmalign + masking
  ├── 04_build_tree.sbatch        # IQ-TREE
  ├── 05_prune_tree.sbatch        # Treemmer pruning
  ├── submit_pipeline.sh          # Master job submitter
  ├── mask_alignment.py           # Gappy-column masking
  └── extract_taxonomy.py         # Taxonomy extraction
models/
  └── RF00177.cm                  # SSU rRNA covariance model (Infernal)
results/
  ├── iqtree/                     # IQ-TREE outputs
  └── treemmer/                   # Treemmer outputs
logs/
  ├── pipeline.log                # Master pipeline log
  └── silva-*.out / silva-*.err   # SLURM job logs
```

## Requirements

**System:**
- Agate supercomputer (UMN HPC) or similar HPC cluster
- SLURM job scheduler
- ~2.2 GB disk space for SILVA input
- ~500 GB disk space for outputs

**Modules (available on Agate):**
```bash
module load vsearch
module load infernal
module load iqtree
module load parallel
module load python3
```

**Python packages** (install via `pip install -r requirements.txt`):
- scikit-bio ≥0.5.7 — FASTA/alignment parsing
- numpy ≥1.19.0 — array operations for masking
- treemmer ≥0.3 — tree pruning with PD maximization

## Notes

### RF00177.cm Scope
The RF00177 covariance model is trained on **bacterial** SSU rRNA sequences. The SILVA NR99 database contains bacteria, archaea, and eukaryotes. The alignment may be suboptimal for archaea and eukaryotic sequences, but they will still align. For a truly comprehensive tree including all domains, consider running separate pipelines with RF01959 (archaea) and RF01960 (eukaryotes).

### Midpoint Rooting
The tree output from IQ-TREE is unrooted. If you need a rooted tree, use scikit-bio:
```bash
python3 -c "
import skbio
tree = skbio.tree.TreeNode.read('results/iqtree/silva_138.2_nr99.treefile', format='newick')
tree_rooted = tree.root_at_midpoint()
with open('results/iqtree/silva_138.2_nr99_rooted.nwk', 'w') as f:
    f.write(str(tree_rooted) + ';')
"
```

## References

- **SILVA 138.2:** https://www.arb-silva.de/
- **VSEARCH:** Rognes et al. (2016). PeerJ 4:e2584
- **cmalign (Infernal):** Nawrocki & Eddy (2013). PLoS Comput Biol 9:e1003213
- **IQ-TREE:** Minh et al. (2020). Mol Biol Evol 37(9):2461–2474
- **Ultrafast bootstrap:** Hoang et al. (2018). Mol Biol Evol 35(2):518–522
- **Treemmer:** Steel & Mooers (2010). Syst Biol 59(6):689–705
- **Ben Kaehler's approach:** https://gist.github.com/BenKaehler/d9291d59bce5cd3d2a90c73b822b3a21

## License

This pipeline code is provided as-is for research purposes. SILVA 138.2 sequences are subject to the SILVA database license terms.
