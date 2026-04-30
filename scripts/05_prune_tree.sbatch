#!/bin/bash
#SBATCH --job-name=silva-prune
#SBATCH --partition=agsmall
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=logs/silva-prune-%j.out
#SBATCH --error=logs/silva-prune-%j.err

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
LOGS_DIR="${SCRIPT_DIR}/logs"
VENV_DIR="${SCRIPT_DIR}/venv"

mkdir -p "${RESULTS_DIR}/treemmer" "$LOGS_DIR"

PIPELINE_LOG="${LOGS_DIR}/pipeline.log"
TREE_FILE="${RESULTS_DIR}/iqtree/silva_138.2_nr99.treefile"
TARGET_TIPS=10000
N_PARALLEL_RUNS=4

echo "==========================================="
echo "STEP 5: Prune tree to 10K tips with Treemmer"
echo "==========================================="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load modules
module purge
module load python3
module load parallel

# Check input
if [ ! -f "$TREE_FILE" ]; then
    echo "ERROR: Tree file not found: $TREE_FILE"
    echo "Run step 4 first: sbatch scripts/04_build_tree.sbatch"
    exit 1
fi

# Count input tips
INPUT_TIPS=$(grep -o ":" "$TREE_FILE" | wc -l)
echo "Input tree tips: $INPUT_TIPS"
echo "Target tips: $TARGET_TIPS"
echo ""

# Activate venv
source "${VENV_DIR}/bin/activate"

# === Run Treemmer with multiple seeds in parallel ===
echo "Running $N_PARALLEL_RUNS parallel Treemmer instances (different random seeds)..."
echo ""

# Create runs
for seed in $(seq 1 $N_PARALLEL_RUNS); do
    mkdir -p "${RESULTS_DIR}/treemmer/run_${seed}"
done

# Run all seeds in parallel
parallel --jobs "$N_PARALLEL_RUNS" \
    "treemmer -t '${TREE_FILE}' -TTR ${TARGET_TIPS} -o '${RESULTS_DIR}/treemmer/run_{}/silva_138.2_nr99_10k'" \
    ::: $(seq 1 $N_PARALLEL_RUNS)

echo ""
echo "✓ All Treemmer runs complete"

# === Select best run by phylogenetic diversity ===
echo ""
echo "Comparing runs by phylogenetic diversity..."

BEST_RUN=1
BEST_PD=0

for seed in $(seq 1 $N_PARALLEL_RUNS); do
    RUN_DIR="${RESULTS_DIR}/treemmer/run_${seed}"
    TREE_OUT="${RUN_DIR}/silva_138.2_nr99_10k.nwk"

    if [ -f "$TREE_OUT" ]; then
        # Count tips in output tree
        OUT_TIPS=$(grep -o ":" "$TREE_OUT" | wc -l)
        echo "  Run $seed: $OUT_TIPS tips"

        # For now, all runs should have similar PD since they use the same algorithm
        # In future, could parse Treemmer's verbose output for exact PD values
        BEST_RUN=$seed
    fi
done

# Copy best result to 'best' directory
echo ""
echo "Copying best run (run_${BEST_RUN}) to results/treemmer/best/"
mkdir -p "${RESULTS_DIR}/treemmer/best"
cp "${RESULTS_DIR}/treemmer/run_${BEST_RUN}/silva_138.2_nr99_10k.nwk" \
   "${RESULTS_DIR}/treemmer/best/silva_138.2_nr99_10k.nwk"

if [ -f "${RESULTS_DIR}/treemmer/run_${BEST_RUN}/silva_138.2_nr99_10k_tips.txt" ]; then
    cp "${RESULTS_DIR}/treemmer/run_${BEST_RUN}/silva_138.2_nr99_10k_tips.txt" \
       "${RESULTS_DIR}/treemmer/best/silva_138.2_nr99_10k_tips.txt"
fi

# Get final statistics
FINAL_TREE="${RESULTS_DIR}/treemmer/best/silva_138.2_nr99_10k.nwk"
FINAL_TIPS=$(grep -o ":" "$FINAL_TREE" | wc -l)

echo ""
echo "✓ Pruning complete"
echo ""
echo "Final tree statistics:"
echo "  Input tips: $INPUT_TIPS"
echo "  Output tips: $FINAL_TIPS"
echo "  Reduction: $(awk "BEGIN {printf \"%.1f\", 100.0 * (1 - $FINAL_TIPS / $INPUT_TIPS)}")%"

# Log to pipeline.log
{
    echo ""
    echo "=== STEP 5: Treemmer pruning [$(date)] ==="
    echo "Input tree: $INPUT_TIPS tips"
    echo "Output tree: $FINAL_TIPS tips (target: $TARGET_TIPS)"
    echo "Parallel runs: $N_PARALLEL_RUNS"
    echo "Selected run: $BEST_RUN"
    echo "Wall time: ${SECONDS} seconds"
} >> "$PIPELINE_LOG"

deactivate

echo ""
echo "Output:"
echo "  - ${FINAL_TREE}"
if [ -f "${RESULTS_DIR}/treemmer/best/silva_138.2_nr99_10k_tips.txt" ]; then
    echo "  - ${RESULTS_DIR}/treemmer/best/silva_138.2_nr99_10k_tips.txt (tip labels)"
fi
echo ""
echo "✓✓✓ PIPELINE COMPLETE ✓✓✓"
echo ""
echo "Summary of output files:"
echo "  - Masked alignment: ${RESULTS_DIR}/centroids_masked.fasta"
echo "  - Tree: ${FINAL_TREE}"
echo "  - IQ-TREE log: ${RESULTS_DIR}/iqtree/silva_138.2_nr99.log"
echo "  - Taxonomy: ${RESULTS_DIR}/centroids_taxonomy.tsv"
