#!/bin/bash
#SBATCH --job-name=silva-tree
#SBATCH --partition=agsmall
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=logs/silva-tree-%j.out
#SBATCH --error=logs/silva-tree-%j.err
set -e

# Explicitly change to repo directory
cd /users/9/falb0011/silva_tree_build/SILVA_138_iqTree/SILVA_138.2_PhylogeneticTree

SCRIPT_DIR="$(pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
LOGS_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${RESULTS_DIR}/fasttree" "$LOGS_DIR"

PIPELINE_LOG="${LOGS_DIR}/pipeline.log"
MASKED_FASTA="${RESULTS_DIR}/centroids_masked.fasta"
TREE_NWK="${RESULTS_DIR}/fasttree/silva_138.2_nr99.nwk"
TREE_LOG="${RESULTS_DIR}/fasttree/silva_138.2_nr99.log"

echo "========================================="
echo "STEP 4: Build tree with FastTree2"
echo "========================================="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load FastTree (multi-threaded OpenMP build)
module purge
module load fasttree

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

echo "Note: module is single-threaded SSE3 build (no OpenMP)"
echo ""

# Run FastTree with GTR+gamma model on nucleotides
# -gtr:     generalized time-reversible model (appropriate for 16S rRNA)
# -gamma:   rescale branch lengths to optimize Gamma20 likelihood after CAT
# -nt:      nucleotide alignment (not amino acid)
# -fastest: faster NJ phase + reduced memory (recommended for >50,000 sequences)
# -log:     detailed log file
echo "Running FastTree (GTR+gamma, nucleotide, fastest mode)..."
fasttree \
    -gtr \
    -gamma \
    -nt \
    -fastest \
    -log "$TREE_LOG" \
    "$MASKED_FASTA" \
    > "$TREE_NWK"

# Check for successful completion
if [ ! -f "$TREE_NWK" ]; then
    echo "ERROR: Tree file not found after FastTree run"
    exit 1
fi

TREE_SIZE=$(wc -c < "$TREE_NWK")
TREE_TIPS=$(grep -o ":" "$TREE_NWK" | wc -l)
echo ""
echo "Tree building complete."
echo "  Output:    $TREE_NWK"
echo "  File size: $TREE_SIZE bytes"
echo "  Tips:      $TREE_TIPS"

# Log to pipeline.log
{
    echo ""
    echo "=== STEP 4: FastTree2 [$(date)] ==="
    echo "Tool:   FastTreeMP (GTR+gamma, nucleotide)"
    echo "Input:  $SEQ_COUNT sequences x ${ALN_WIDTH} columns"
    echo "Output: $TREE_TIPS tips"
    echo "Wall time: ${SECONDS} seconds"
} >> "$PIPELINE_LOG"

echo ""
echo "✓ Tree building complete"
echo "Outputs:"
echo "  - $TREE_NWK (Newick tree)"
echo "  - $TREE_LOG (FastTree log)"
echo ""
echo "Next step:"
echo "  sbatch scripts/05_prune_tree.sbatch"
