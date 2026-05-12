# SILVA 138.2 NR99 Bacterial Phylogenetic Tree Pipeline

Builds a full-resolution bacterial reference phylogenetic tree from SILVA NR99 138.2, suitable for phylogenetic placement of 16S amplicon data using pplacer-BSCAMPP.

## Overview

The pipeline runs as five modular SLURM jobs on Agate (UMN HPC):

1. **Download** the unaligned SILVA NR99 138.2 FASTA (~510K sequences)
2. **Cluster** at 97% identity (VSEARCH), then filter to bacteria only — yields 167,957 centroids
3. **Align & mask** centroid sequences to RF00177.cm (cmalign), convert, and mask gappy columns
4. **Build tree** with FastTree2 (GTR+gamma, `-fastest`) — full 167,957-tip bacterial tree
5. **Set up SCAMPP** — create a Taxtastic reference package for pplacer-BSCAMPP placement

No tree pruning is performed. The full 167,957-tip tree is used directly with pplacer-BSCAMPP, which handles backbone trees up to 200,000 leaves (Wedell, Cai & Warnow 2023).

### Final Output Files

| File | Description |
|------|-------------|
| `results/centroids_masked.fasta` | Masked alignment used for tree building (~1533 columns) |
| `results/fasttree/silva_138.2_nr99.nwk` | Full bacterial tree, 167,957 tips (Newick) |
| `results/fasttree/silva_138.2_nr99.log` | FastTree log (GTR+gamma model parameters) |
| `results/centroids_taxonomy.tsv` | Full SILVA lineage per centroid (seqid + semicolon-delimited lineage) |
| `results/scampp/silva_138.2_nr99.refpkg/` | Taxtastic reference package for pplacer |
| `logs/pipeline.log` | Per-step statistics summary |

## Setup (once, interactive on Agate login node)

### 1. Download SILVA

```bash
cd /path/to/SILVA_138.2_NR99_fastTree
bash scripts/01_download_silva.sh
```

Downloads `SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz` (~2.2 GB). Safe to re-run (checks for existing file).

### 2. Create Python virtual environment

```bash
module load python3          # loads Python 3.10+
python3 -m venv venv
source venv/bin/activate
pip install -U pip
pip install -r requirements.txt
deactivate
```

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
3. **mask_alignment.py** — removes columns with ≥99.56% gaps (numpy, fast)
4. **extract_taxonomy.py** — parses SILVA FASTA headers to extract full lineage

**Output:**
- `results/centroids_aligned.fasta` — aligned FASTA (1533 columns)
- `results/centroids_masked.fasta` — masked alignment (~1533 columns, gappy columns removed)
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
- `results/fasttree/silva_138.2_nr99.log` — FastTree log with GTR+gamma parameters

### Step 5: SCAMPP Setup (Reference Package)

**File:** [scripts/05_setup_scampp.sbatch](scripts/05_setup_scampp.sbatch)  
**Resources:** agsmall, 1 core, 32 GB, 2 hours  
**Module:** `python3`

No pruning is performed. The full 167,957-tip tree is within pplacer-BSCAMPP's tested range (up to 200,000 leaves; Wedell et al. 2023).

1. Installs `bscampp` and `taxtastic` into the venv
2. Creates a Taxtastic reference package from the FastTree tree and its log file (which provides GTR+gamma model parameters for pplacer)
3. Writes a per-study placement script: `scripts/run_placement.sh`

```bash
taxit create \
    -l "silva_138.2_nr99_bacteria" \
    -P results/scampp/silva_138.2_nr99.refpkg \
    --aln-fasta results/centroids_masked.fasta \
    --tree results/fasttree/silva_138.2_nr99.nwk \
    --tree-stats results/fasttree/silva_138.2_nr99.log
```

**Output:**
- `results/scampp/silva_138.2_nr99.refpkg/` — Taxtastic reference package
- `scripts/run_placement.sh` — template script for per-study placement

## Placing Query Sequences

For each study with 16S ASVs or OTUs:

```bash
bash scripts/run_placement.sh /path/to/query_reads.fasta results/scampp/my_study 8
```

The script does two things:
1. **Aligns queries** to RF00177.cm with cmalign (`--matchonly --dnaout`)
2. **Places queries** into the reference tree using pplacer-BSCAMPP (`-b 2000`)

Output: `results/scampp/my_study.jplace` — standard placement file readable by guppy, iTOL, etc.

**Summarize placements:**
```bash
guppy fat results/scampp/my_study.jplace
```

**Note on pplacer binary:** pplacer must be available on PATH before running placements. Download from [https://github.com/matsen/pplacer/releases](https://github.com/matsen/pplacer/releases) and place in `bin/`. pplacer standalone (without SCAMPP) cannot handle trees >78K tips; always use `run_bscampp.py` (installed with `pip install bscampp`).

## Repository Structure

```
scripts/
  01_download_silva.sh        Download SILVA NR99 138.2 FASTA
  02_cluster.sbatch           VSEARCH 97% clustering + bacteria filter
  03_align_and_mask.sbatch    cmalign + pfam_to_fasta + masking + taxonomy
  04_build_tree.sbatch        FastTree2 GTR+gamma
  05_setup_scampp.sbatch      Install BSCAMPP + create Taxtastic refpkg
  run_placement.sh            Per-study placement (generated by step 5)
  submit_pipeline.sh          Submit all jobs with dependency chaining
  filter_bacteria.py          Filter centroids to bacteria only
  pfam_to_fasta.py            Convert Pfam/Stockholm to FASTA
  mask_alignment.py           Remove gappy alignment columns
  extract_taxonomy.py         Extract SILVA taxonomy from FASTA headers
models/
  RF00177.cm                  Bacteria SSU rRNA covariance model (Infernal)
results/
  centroids_97_bacteria.fasta Bacteria centroids
  centroids_masked.fasta      Masked alignment (tree input)
  fasttree/                   FastTree outputs (tree + log)
  scampp/                     Reference package + per-study placements
  centroids_taxonomy.tsv      Taxonomy per tip
logs/
  pipeline.log                Master pipeline log
  silva-*.out / *.err         SLURM job logs
```

## Troubleshooting

**Step 2 (clustering) fails:**  
Check `logs/silva-cluster-*.err`. Verify SILVA FASTA exists: `ls -lh SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz`

**Step 3 (cmalign) fails with out-of-memory:**  
Increase `--mxsize` in the cmalign call (try `16000`) or switch to `aglarge` partition.

**Step 4 (FastTree) produces empty tree:**  
Check `logs/silva-tree-*.err`. FastTree writes the tree to stdout; ensure the redirect is intact. Verify masked FASTA is not empty: `grep -c "^>" results/centroids_masked.fasta`

**Step 5 (taxit) fails:**  
The FastTree log (`silva_138.2_nr99.log`) must exist and contain GTR rate parameters. Confirm step 4 completed: `ls -lh results/fasttree/`

**run_placement.sh fails at cmalign:**  
Ensure `module load infernal` succeeds in your environment. RF00177.cm must be at `models/RF00177.cm`.

## Requirements

**HPC modules (available on Agate):**
```
vsearch     Step 2: clustering
infernal    Step 3: cmalign
fasttree    Step 4: tree building
python3     Steps 3, 5, and setup
```

**Python packages** (installed into `venv/` via `pip install -r requirements.txt`):
- `numpy` — array operations for alignment masking
- `bscampp` — installed in step 5 via pip
- `taxtastic` — installed in step 5 via pip

## References

- **SILVA 138.2:** Quast et al. (2013). Nucleic Acids Res 41:D590–D596. https://www.arb-silva.de/
- **VSEARCH:** Rognes et al. (2016). PeerJ 4:e2584
- **Infernal/cmalign:** Nawrocki & Eddy (2013). PLoS Comput Biol 9:e1003213
- **FastTree2:** Price, Dehal & Arkin (2010). PLoS ONE 5(3):e9490
- **pplacer-BSCAMPP:** Wedell, Cai & Warnow (2023). IEEE/ACM TCBB 20(2):1417–1430
- **RF00177:** Rfam database — bacteria SSU rRNA covariance model
