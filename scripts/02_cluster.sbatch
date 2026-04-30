#!/bin/bash
#SBATCH --job-name=silva-cluster
#SBATCH --partition=agsmall
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=logs/silva-cluster-%j.out
#SBATCH --error=logs/silva-cluster-%j.err

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SILVA_FILE="${SCRIPT_DIR}/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
RESULTS_DIR="${SCRIPT_DIR}/results"
LOGS_DIR="${SCRIPT_DIR}/logs"

mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

PIPELINE_LOG="${LOGS_DIR}/pipeline.log"

echo "=========================================="
echo "STEP 2: Cluster at 97% similarity (VSEARCH)"
echo "=========================================="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load vsearch
module purge
module load vsearch

# Check input file
if [ ! -f "$SILVA_FILE" ]; then
    echo "ERROR: Input file not found: $SILVA_FILE"
    echo "Run: bash scripts/01_download_silva.sh"
    exit 1
fi

# Count input sequences
echo "Counting input sequences..."
if [[ "$SILVA_FILE" == *.gz ]]; then
    INPUT_COUNT=$(zcat "$SILVA_FILE" | grep -c "^>")
else
    INPUT_COUNT=$(grep -c "^>" "$SILVA_FILE")
fi
echo "Input sequences: $INPUT_COUNT"

# Run VSEARCH clustering
echo ""
echo "Running VSEARCH 97% clustering..."
echo "Command: vsearch --cluster_fast <input> --id 0.97 --threads $SLURM_CPUS_PER_TASK --centroids <out>"
echo ""

vsearch --cluster_fast "$SILVA_FILE" \
    --id 0.97 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --centroids "${RESULTS_DIR}/centroids_97.fasta" \
    --uc "${RESULTS_DIR}/clusters_97.uc" \
    --fasta_width 0

# Count output centroids
echo ""
CENTROID_COUNT=$(grep -c "^>" "${RESULTS_DIR}/centroids_97.fasta")
echo "Output centroids: $CENTROID_COUNT"

# Compression stats
CLUSTER_RATIO=$(awk "BEGIN {printf \"%.1f\", 100.0 * $CENTROID_COUNT / $INPUT_COUNT}")
echo "Cluster compression: $CLUSTER_RATIO%"

# Log to pipeline.log
{
    echo ""
    echo "=== STEP 2: Clustering [$(date)] ==="
    echo "Input sequences: $INPUT_COUNT"
    echo "Output centroids: $CENTROID_COUNT"
    echo "Cluster compression: ${CLUSTER_RATIO}%"
    echo "Identity threshold: 97%"
    echo "Wall time: ${SECONDS} seconds"
} >> "$PIPELINE_LOG"

echo ""
echo "✓ Clustering complete"
echo "Outputs:"
echo "  - ${RESULTS_DIR}/centroids_97.fasta"
echo "  - ${RESULTS_DIR}/clusters_97.uc"
echo ""
echo "Next step:"
echo "  sbatch scripts/03_align_and_mask.sbatch"
