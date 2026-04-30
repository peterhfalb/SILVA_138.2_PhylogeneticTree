#!/bin/bash
# Submit all pipeline jobs with SLURM dependencies
#
# Usage: bash scripts/submit_pipeline.sh
#
# This will:
# 1. Submit step 2 (clustering)
# 2. Submit step 3 (alignment & masking) -- waits for step 2
# 3. Submit step 4 (tree building) -- waits for step 3
# 4. Submit step 5 (pruning) -- waits for step 4
#
# All jobs are chained with --dependency=afterok so they only run if previous succeeds

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOGS_DIR="${SCRIPT_DIR}/logs"

mkdir -p "$LOGS_DIR"

echo "=========================================="
echo "Submitting SILVA phylogenetic tree pipeline"
echo "=========================================="
echo ""

# Initialize pipeline log
PIPELINE_LOG="${LOGS_DIR}/pipeline.log"
{
    echo "========================================"
    echo "SILVA 138.2 NR99 Phylogenetic Tree Pipeline"
    echo "Started: $(date)"
    echo "========================================"
} > "$PIPELINE_LOG"

# === STEP 2: Clustering ===
echo "Submitting STEP 2: Clustering (VSEARCH 97%)..."
JOB_2=$(sbatch -p scripts/02_cluster.sbatch | awk '{print $NF}')
echo "  Job ID: $JOB_2"

# === STEP 3: Alignment & Masking ===
echo ""
echo "Submitting STEP 3: Alignment & Masking (cmalign)..."
echo "  Dependency: afterok:$JOB_2"
JOB_3=$(sbatch --dependency=afterok:$JOB_2 scripts/03_align_and_mask.sbatch | awk '{print $NF}')
echo "  Job ID: $JOB_3"

# === STEP 4: Tree Building ===
echo ""
echo "Submitting STEP 4: Tree Building (IQ-TREE)..."
echo "  Dependency: afterok:$JOB_3"
JOB_4=$(sbatch --dependency=afterok:$JOB_3 scripts/04_build_tree.sbatch | awk '{print $NF}')
echo "  Job ID: $JOB_4"

# === STEP 5: Pruning ===
echo ""
echo "Submitting STEP 5: Pruning to 10K tips (Treemmer)..."
echo "  Dependency: afterok:$JOB_4"
JOB_5=$(sbatch --dependency=afterok:$JOB_4 scripts/05_prune_tree.sbatch | awk '{print $NF}')
echo "  Job ID: $JOB_5"

# === Summary ===
echo ""
echo "=========================================="
echo "Pipeline submitted successfully"
echo "=========================================="
echo ""
echo "Dependency chain:"
echo "  STEP 2 (Clustering):       Job $JOB_2"
echo "  STEP 3 (Alignment):        Job $JOB_3  ← afterok:$JOB_2"
echo "  STEP 4 (IQ-TREE):          Job $JOB_4  ← afterok:$JOB_3"
echo "  STEP 5 (Treemmer):         Job $JOB_5  ← afterok:$JOB_4"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  slog $JOB_5  # Shows final job status when complete"
echo ""
echo "Logs:"
echo "  Pipeline log: ${PIPELINE_LOG}"
echo "  SLURM logs: ${LOGS_DIR}/silva-*.out"
echo ""
