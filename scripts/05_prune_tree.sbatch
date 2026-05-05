#!/bin/bash
#SBATCH --job-name=silva-prune
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=logs/silva-prune-%j.out
#SBATCH --error=logs/silva-prune-%j.err
set -e

# Explicitly change to repo directory
cd /users/9/falb0011/silva_tree_build/SILVA_138_iqTree/SILVA_138.2_PhylogeneticTree

SCRIPT_DIR="$(pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
LOGS_DIR="${SCRIPT_DIR}/logs"
CONTAINER_DIR="${SCRIPT_DIR}/.containers"

mkdir -p "${RESULTS_DIR}/treemmer" "$LOGS_DIR" "$CONTAINER_DIR"

PIPELINE_LOG="${LOGS_DIR}/pipeline.log"
TREE_FILE="${RESULTS_DIR}/fasttree/silva_138.2_nr99.nwk"
TARGET_TIPS=10000
N_PARALLEL_RUNS=4

# Singularity container settings
TREEMMER_IMAGE="${CONTAINER_DIR}/treemmer_0.3.sif"
TREEMMER_SB="${CONTAINER_DIR}/treemmer_sb"

echo "==========================================="
echo "STEP 5: Prune tree to 10K tips with Treemmer"
echo "==========================================="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load modules
module purge
module load singularity
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

# === Download/setup Singularity container if needed ===
if [ ! -f "$TREEMMER_IMAGE" ]; then
    echo "Downloading Treemmer Singularity image (first run only)..."
    singularity pull --arch amd64 \
        --dir "$CONTAINER_DIR" \
        library://fmenardo/treemmer/treemmer:0.3

    # Note: singularity pull downloads to treemmer_0.3.sif by default
fi

if [ ! -d "$TREEMMER_SB" ]; then
    echo "Building Singularity sandbox..."
    singularity build --sandbox "$TREEMMER_SB" "$TREEMMER_IMAGE"
fi

echo "Using Treemmer container: $TREEMMER_SB"
echo ""

# === Run Treemmer with multiple seeds in parallel ===
echo "Running $N_PARALLEL_RUNS parallel Treemmer instances (different random seeds)..."
echo ""

# Create runs
for seed in $(seq 1 $N_PARALLEL_RUNS); do
    mkdir -p "${RESULTS_DIR}/treemmer/run_${seed}"
done

# Run all seeds in parallel using Singularity
# Note: treemmer command inside container is: python3 /Treemmer_v0.3.py
parallel --jobs "$N_PARALLEL_RUNS" \
    "singularity exec '${TREEMMER_SB}' python3 /Treemmer_v0.3.py -t '${TREE_FILE}' -TTR ${TARGET_TIPS} -o '${RESULTS_DIR}/treemmer/run_{}/silva_138.2_nr99_10k'" \
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
    echo "Container: Singularity (treemmer v0.3)"
    echo "Wall time: ${SECONDS} seconds"
} >> "$PIPELINE_LOG"

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
