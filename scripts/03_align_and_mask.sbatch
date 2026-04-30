#!/bin/bash
#SBATCH --job-name=silva-align
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=logs/silva-align-%j.out
#SBATCH --error=logs/silva-align-%j.err
set -e

# Explicitly change to repo directory
cd /users/9/falb0011/silva_tree_build/SILVA_138_iqTree/SILVA_138.2_PhylogeneticTree

SCRIPT_DIR="$(pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
LOGS_DIR="${SCRIPT_DIR}/logs"
MODELS_DIR="${SCRIPT_DIR}/models"
VENV_DIR="${SCRIPT_DIR}/venv"

mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

PIPELINE_LOG="${LOGS_DIR}/pipeline.log"
CENTROID_FASTA="${RESULTS_DIR}/centroids_97_bacteria.fasta"
ALIGNED_FASTA="${RESULTS_DIR}/centroids_aligned.fasta"
MASKED_FASTA="${RESULTS_DIR}/centroids_masked.fasta"
TAXONOMY_TSV="${RESULTS_DIR}/centroids_taxonomy.tsv"

echo "==========================================="
echo "STEP 3: Align & Mask (cmalign + masking)"
echo "==========================================="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load modules
module purge
module load infernal/1.1.4
module load python3

# Check input
if [ ! -f "$CENTROID_FASTA" ]; then
    echo "ERROR: Centroid FASTA not found: $CENTROID_FASTA"
    echo "Run step 2 first: sbatch scripts/02_cluster.sbatch"
    exit 1
fi

if [ ! -f "${MODELS_DIR}/RF00177.cm" ]; then
    echo "ERROR: Covariance model not found: ${MODELS_DIR}/RF00177.cm"
    exit 1
fi

# === Alignment with cmalign ===
echo "Step 3a: Aligning to RF00177.cm..."
CENTROID_COUNT=$(grep -c "^>" "$CENTROID_FASTA")
echo "  Input sequences: $CENTROID_COUNT"

# Use Pfam format (can handle large alignments) then convert to FASTA
ALIGNED_PFAM="${RESULTS_DIR}/centroids_aligned.pfam"

cmalign --cpu "$SLURM_CPUS_PER_TASK" \
    --mxsize 8000 \
    --outformat Pfam \
    "${MODELS_DIR}/RF00177.cm" \
    "$CENTROID_FASTA" \
    > "$ALIGNED_PFAM"

# Convert Pfam format to FASTA using esl-reformat
esl-reformat fasta "$ALIGNED_PFAM" > "$ALIGNED_FASTA"
rm "$ALIGNED_PFAM"

ALIGNED_COUNT=$(grep -c "^>" "$ALIGNED_FASTA")
echo "  Aligned sequences: $ALIGNED_COUNT"

# Get alignment dimensions
ALIGNMENT_WIDTH=$(grep -v "^>" "$ALIGNED_FASTA" | head -1 | wc -c)
echo "  Alignment width: $ALIGNMENT_WIDTH columns"

# === Masking ===
echo ""
echo "Step 3b: Masking gappy columns..."
source "${VENV_DIR}/bin/activate"
python3 "${SCRIPT_DIR}/scripts/mask_alignment.py" \
    "$ALIGNED_FASTA" \
    "$MASKED_FASTA" \
    --gap-threshold 99.56

MASKED_COUNT=$(grep -c "^>" "$MASKED_FASTA")
MASKED_WIDTH=$(grep -v "^>" "$MASKED_FASTA" | head -1 | wc -c)
echo "  Masked alignment: $MASKED_COUNT sequences × $MASKED_WIDTH columns"

# === Taxonomy extraction ===
echo ""
echo "Step 3c: Extracting taxonomy..."
python3 "${SCRIPT_DIR}/scripts/extract_taxonomy.py" \
    "$CENTROID_FASTA" \
    "$TAXONOMY_TSV"

TAXONOMY_LINES=$(wc -l < "$TAXONOMY_TSV")
echo "  Taxonomy entries: $((TAXONOMY_LINES - 1)) (excluding header)"

deactivate

# === Log results ===
echo ""
{
    echo ""
    echo "=== STEP 3: Align & Mask [$(date)] ==="
    echo "Alignment:"
    echo "  Input: $CENTROID_COUNT sequences"
    echo "  Output: $ALIGNED_COUNT sequences × ~${ALIGNMENT_WIDTH} columns"
    echo "Masking:"
    echo "  Output: $MASKED_COUNT sequences × $MASKED_WIDTH columns"
    echo "  Threshold: 99.56% gaps"
    echo "Taxonomy:"
    echo "  Entries: $((TAXONOMY_LINES - 1))"
    echo "Wall time: ${SECONDS} seconds"
} >> "$PIPELINE_LOG"

echo ""
echo "✓ Align & Mask complete"
echo "Outputs:"
echo "  - $ALIGNED_FASTA"
echo "  - $MASKED_FASTA"
echo "  - $TAXONOMY_TSV"
echo ""
echo "Next step:"
echo "  sbatch scripts/04_build_tree.sbatch"
