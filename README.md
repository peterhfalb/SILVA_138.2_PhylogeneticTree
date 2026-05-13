# SILVA 138.2 NR99 Bacterial Phylogenetic Tree Pipeline

Builds a full-resolution bacterial reference phylogenetic tree from SILVA NR99 138.2, suitable for phylogenetic placement of 16S amplicon data using pplacer-BSCAMPP.

## Overview

The pipeline runs as five modular SLURM jobs on Agate (UMN HPC):

1. **Download** the unaligned SILVA NR99 138.2 FASTA (~510K sequences)
2. **Cluster** at 97% identity (VSEARCH), then filter to bacteria only — yields 167,957 centroids
3. **Align & mask** centroid sequences to RF00177.cm (cmalign), convert, and mask gappy columns; saves the column mask for use during placement
4. **Build tree** with FastTree2 (GTR+gamma, `-fastest`) — full 167,957-tip bacterial tree
5. **Set up BSCAMPP** — install bscampp and write the per-study placement script

No tree pruning is performed. The full 167,957-tip tree is used directly with pplacer-BSCAMPP, which handles backbone trees up to 200,000 leaves (Wedell, Cai & Warnow 2023).

### Final Output Files

| File | Description |
|------|-------------|
| `results/centroids_masked.fasta` | Masked reference alignment used to build the tree |
| `results/centroids_column_mask.npy` | Boolean column mask — which columns were retained after masking |
| `results/fasttree/silva_138.2_nr99.nwk` | Full bacterial tree, 167,957 tips (Newick) |
| `results/fasttree/silva_138.2_nr99.log` | FastTree log (GTR+gamma model parameters, read directly by bscampp) |
| `results/centroids_taxonomy.tsv` | Full SILVA lineage per centroid (seqid + semicolon-delimited lineage) |
| `logs/pipeline.log` | Per-step statistics summary |

## Setup (once, interactive on Agate login node)

### 1. Create Python virtual environment

```bash
module load python3          # loads Python 3.10+
python3 -m venv venv
source venv/bin/activate
pip install -U pip
pip install -r requirements.txt
deactivate
```

### 2. Download SILVA

```bash
cd /path/to/SILVA_138.2_NR99_fastTree
bash scripts/01_download_silva.sh
```

Downloads `SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz` (~2.2 GB). Safe to re-run (skips files that already exist).

## Running the Pipeline

### Option A: Submit all jobs with automatic dependencies

```bash
bash scripts/submit_pipeline.sh
```

Submits jobs 2–5 chained with `--dependency=afterok`, so each step only starts if the previous succeeded.

### Option B: Submit steps individually

```bash
sbatch scripts/02_cluster.sbatch
sbatch scripts/03_align_and_mask.sbatch   # after step 2
sbatch scripts/04_build_tree.sbatch       # after step 3
sbatch scripts/05_setup_scampp.sbatch     # after step 4
```

Monitor: `squeue -u $USER`

## Pipeline Details

### Step 2: Clustering & Bacteria Filtering

**File:** [scripts/02_cluster.sbatch](scripts/02_cluster.sbatch)  
**Resources:** agsmall, 32 cores, 128 GB, 4 hours  
**Modules:** `vsearch`, `python3`

1. VSEARCH 97% clustering: `--cluster_fast --id 0.97 --notrunclabels`
2. `filter_bacteria.py` — keeps only bacterial centroids using SILVA taxonomy headers

**Output:**
- `results/centroids_97.fasta` — all centroids (206,289 total)
- `results/centroids_97_bacteria.fasta` — bacteria-only centroids (167,957)
- `results/clusters_97.uc` — cluster assignments

### Step 3: Alignment & Masking

**File:** [scripts/03_align_and_mask.sbatch](scripts/03_align_and_mask.sbatch)  
**Resources:** agsmall, 16 cores, 128 GB, 12 hours  
**Modules:** `infernal`, `python3`

1. **cmalign** — aligns bacteria centroids to RF00177.cm (bacteria SSU rRNA model)
   - `--matchonly` — outputs only the 1533 consensus match columns (drops insert columns)
   - `--dnaout` — forces DNA output (T not U)
   - `--outformat Pfam` — multi-block Stockholm format (required for large alignments)
   - `--mxsize 12000 --maxtau 0.2` — memory and accuracy tuning for 167K sequences
2. **pfam_to_fasta.py** — converts Pfam/Stockholm multi-block output to FASTA
3. **mask_alignment.py** — removes columns with ≥99.56% gaps; saves the retained column indices to `centroids_column_mask.npy`
4. **extract_taxonomy.py** — parses SILVA FASTA headers to extract full lineage

**Output:**
- `results/centroids_aligned.fasta` — aligned FASTA (1533 columns)
- `results/centroids_masked.fasta` — masked alignment (gappy columns removed)
- `results/centroids_column_mask.npy` — boolean array of retained column positions (applied to query alignments during placement)
- `results/centroids_taxonomy.tsv` — `seqid\ttaxonomy` (full semicolon-delimited SILVA lineage)

### Step 4: Tree Building (FastTree2)

**File:** [scripts/04_build_tree.sbatch](scripts/04_build_tree.sbatch)  
**Resources:** agsmall, 1 core, 64 GB, 36–48 hours  
**Module:** `fasttree` (single-threaded SSE3 build; command: `fasttree`)

```bash
fasttree -gtr -gamma -nt -fastest -log "$TREE_LOG" "$MASKED_FASTA" > "$TREE_NWK"
```

- `-gtr` — generalized time-reversible model (standard for 16S rRNA)
- `-gamma` — rescales branch lengths to optimize Gamma20 likelihood after CAT approximation
- `-nt` — nucleotide alignment
- `-fastest` — required for >50,000 sequences (faster NJ phase, reduced memory)

**Note:** IQ-TREE was evaluated but found infeasible at 167K sequences — thread benchmarking showed 39% efficiency at 2 threads, and O(N²) scaling made completion impossible within wall time limits. FastTree2 with `-fastest` is the standard approach for large 16S reference trees (used by SILVA itself).

**Output:**
- `results/fasttree/silva_138.2_nr99.nwk` — full 167,957-tip tree (Newick)
- `results/fasttree/silva_138.2_nr99.log` — FastTree log with GTR+gamma parameters (used directly by bscampp as `-i`)

### Step 5: BSCAMPP Setup

**File:** [scripts/05_setup_scampp.sbatch](scripts/05_setup_scampp.sbatch)  
**Resources:** agsmall, 1 core, 8 GB, 30 minutes  
**Module:** `python3`

Installs bscampp (`pip install bscampp`) and writes `scripts/run_placement.sh` — a template for per-study query placement. No reference package is built; bscampp takes the tree, alignment, and FastTree log as direct file arguments.

## Placing Query Sequences

For each study with 16S ASVs or OTUs, run the generated placement script:

```bash
bash scripts/run_placement.sh query_reads.fasta results/placements my_study 8
```

The script performs three steps:

1. **Align queries** to RF00177.cm with cmalign (`--matchonly --dnaout`) → 1533 match columns
2. **Apply column mask** (`--apply-mask centroids_column_mask.npy`) — trims query alignment to the exact same columns used to build the reference tree; reference and query **must** have identical columns for pplacer
3. **Place** with pplacer-BSCAMPP:

```bash
run_bscampp.py \
    -i results/fasttree/silva_138.2_nr99.log \
    -t results/fasttree/silva_138.2_nr99.nwk \
    -a results/centroids_masked.fasta \
    -q query_aligned_masked.fasta \
    -d results/placements \
    -o my_study \
    --threads 8
```

Output: `results/placements/my_study.jplace` — standard placement file readable by guppy, iTOL, etc.

## Repository Structure

```
scripts/
  01_download_silva.sh      Download SILVA NR99 138.2 FASTA
  02_cluster.sbatch         VSEARCH 97% clustering + bacteria filter
  03_align_and_mask.sbatch  cmalign + pfam_to_fasta + masking + taxonomy
  04_build_tree.sbatch      FastTree2 GTR+gamma
  05_setup_scampp.sbatch    Install bscampp + write placement script
  run_placement.sh          Per-study placement (generated by step 5)
  submit_pipeline.sh        Submit all jobs with dependency chaining
  filter_bacteria.py        Filter centroids to bacteria only
  pfam_to_fasta.py          Convert Pfam/Stockholm to FASTA
  mask_alignment.py         Remove gappy alignment columns; save/apply column mask
  extract_taxonomy.py       Extract SILVA taxonomy from FASTA headers
models/
  RF00177.cm                Bacteria SSU rRNA covariance model (Infernal)
results/
  centroids_97_bacteria.fasta   Bacteria centroids
  centroids_masked.fasta        Masked alignment (reference for tree + placement)
  centroids_column_mask.npy     Column mask applied to query alignments
  centroids_taxonomy.tsv        Taxonomy per tip
  fasttree/                     FastTree tree + log
  placements/                   Per-study jplace output files
logs/
  pipeline.log              Master pipeline log
  silva-*.out / *.err       SLURM job logs
```

## Troubleshooting

**Step 2 (clustering) fails:**  
Check `logs/silva-cluster-*.err`. Verify SILVA FASTA exists: `ls -lh SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz`

**Step 3 (cmalign) fails with out-of-memory:**  
Increase `--mxsize` in the cmalign call (try `16000`) or switch to `aglarge` partition.

**Step 4 (FastTree) produces empty tree:**  
Check `logs/silva-tree-*.err`. FastTree writes the tree to stdout; ensure the redirect is intact. Verify masked FASTA is not empty: `grep -c "^>" results/centroids_masked.fasta`

**Step 5 (bscampp install) fails:**  
Ensure the venv is activated and Python 3.10+ is loaded. Try `pip install -U bscampp`.

**Placement fails with column mismatch:**  
The query alignment must have the same number of columns as `centroids_masked.fasta`. Ensure `--apply-mask centroids_column_mask.npy` is being used in `run_placement.sh`, not a fresh gap-masking run on the query sequences.

**Placement fails at cmalign:**  
Ensure `module load infernal` succeeds. RF00177.cm must be at `models/RF00177.cm`.

## Requirements

**HPC modules (available on Agate):**
```
vsearch     Step 2: clustering
infernal    Step 3: cmalign + query alignment
fasttree    Step 4: tree building
python3     Steps 3, 5, setup, and placement
```

**Python packages** (installed into `venv/` via `pip install -r requirements.txt` + step 5):
- `scikit-bio` — alignment I/O in masking script
- `numpy` — array operations for column masking
- `bscampp` — installed in step 5 via `pip install bscampp`

## References

- **SILVA 138.2:** Quast et al. (2013). Nucleic Acids Res 41:D590–D596. https://www.arb-silva.de/
- **VSEARCH:** Rognes et al. (2016). PeerJ 4:e2584
- **Infernal/cmalign:** Nawrocki & Eddy (2013). PLoS Comput Biol 9:e1003213
- **FastTree2:** Price, Dehal & Arkin (2010). PLoS ONE 5(3):e9490
- **pplacer-BSCAMPP:** Wedell, Cai & Warnow (2023). IEEE/ACM TCBB 20(2):1417–1430
- **RF00177:** Rfam database — bacteria SSU rRNA covariance model
