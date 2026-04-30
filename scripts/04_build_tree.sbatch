#!/bin/bash
#SBATCH --job-name=silva-tree
#SBATCH --partition=aglarge
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --output=logs/silva-tree-%j.out
#SBATCH --error=logs/silva-tree-%j.err

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
LOGS_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${RESULTS_DIR}/iqtree" "$LOGS_DIR"

PIPELINE_LOG="${LOGS_DIR}/pipeline.log"
MASKED_FASTA="${RESULTS_DIR}/centroids_masked.fasta"
TREE_PREFIX="${RESULTS_DIR}/iqtree/silva_138.2_nr99"

echo "========================================="
echo "STEP 4: Build tree with IQ-TREE"
echo "========================================="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load iqtree
module purge
module load iqtree

# Check input
if [ ! -f "$MASKED_FASTA" ]; then
    echo "ERROR: Masked alignment not found: $MASKED_FASTA"
    echo "Run step 3 first: sbatch scripts/03_align_and_mask.sbatch"
    exit 1
fi

SEQ_COUNT=$(grep -c "^>" "$MASKED_FASTA")
echo "Input sequences: $SEQ_COUNT"
ALN_WIDTH=$(grep -v "^>" "$MASKED_FASTA" | head -1 | wc -c)
echo "Alignment width: $ALN_WIDTH columns"
echo ""

# Run IQ-TREE
echo "Running IQ-TREE2 with ModelFinder + 1000 ultrafast bootstrap replicates..."
echo "Command: iqtree2 -s $MASKED_FASTA -m TEST -bb 1000 -nt AUTO -ntmax 64"
echo ""

iqtree2 \
    -s "$MASKED_FASTA" \
    -m TEST \
    -bb 1000 \
    -nt AUTO \
    -ntmax "$SLURM_CPUS_PER_TASK" \
    -pre "$TREE_PREFIX" \
    -seed 12345

# Check for successful completion
if [ ! -f "${TREE_PREFIX}.treefile" ]; then
    echo "ERROR: Tree file not found after IQ-TREE run"
    exit 1
fi

# Extract key statistics from IQ-TREE output
echo ""
echo "Tree building complete. Key outputs:"
ls -lh "${TREE_PREFIX}".{treefile,log,iqtree}

# Parse model from .iqtree file
echo ""
echo "Model selection results:"
if [ -f "${TREE_PREFIX}.iqtree" ]; then
    grep "Best-fit model" "${TREE_PREFIX}.iqtree" | head -1
fi

TREE_TIPS=$(grep -o ":" "${TREE_PREFIX}.treefile" | wc -l)
echo "Tree tips: $TREE_TIPS"

# Log to pipeline.log
IQTREE_LOG="${TREE_PREFIX}.log"
BEST_MODEL=$(grep "Best-fit model" "${TREE_PREFIX}.iqtree" | head -1 | sed 's/^.*: //' || echo "unknown")

{
    echo ""
    echo "=== STEP 4: IQ-TREE [$(date)] ==="
    echo "Input: $SEQ_COUNT sequences × ${ALN_WIDTH} columns"
    echo "Model selection: -m TEST (ModelFinder)"
    echo "Bootstrap: -bb 1000 (ultrafast)"
    echo "Best model: $BEST_MODEL"
    echo "Output tips: $TREE_TIPS"
    echo "Wall time: ${SECONDS} seconds"
} >> "$PIPELINE_LOG"

echo ""
echo "✓ Tree building complete"
echo "Outputs:"
echo "  - ${TREE_PREFIX}.treefile (ML tree)"
echo "  - ${TREE_PREFIX}.log (IQ-TREE log)"
echo "  - ${TREE_PREFIX}.iqtree (full report)"
echo ""
echo "Next step:"
echo "  sbatch scripts/05_prune_tree.sbatch"
