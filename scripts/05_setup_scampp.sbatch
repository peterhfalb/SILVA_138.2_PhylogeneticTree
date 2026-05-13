#!/bin/bash
#SBATCH --job-name=silva-scampp-setup
#SBATCH --partition=agsmall
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=logs/silva-scampp-setup-%j.out
#SBATCH --error=logs/silva-scampp-setup-%j.err
set -e

# Explicitly change to repo directory
cd /users/9/falb0011/silva_tree_build/SILVA_138_iqTree/SILVA_138.2_PhylogeneticTree

SCRIPT_DIR="$(pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
LOGS_DIR="${SCRIPT_DIR}/logs"
VENV_DIR="${SCRIPT_DIR}/venv"

mkdir -p "$LOGS_DIR"

PIPELINE_LOG="${LOGS_DIR}/pipeline.log"
TREE_FILE="${RESULTS_DIR}/fasttree/silva_138.2_nr99.nwk"
TREE_LOG="${RESULTS_DIR}/fasttree/silva_138.2_nr99.log"
ALN_FILE="${RESULTS_DIR}/centroids_masked.fasta"
COLUMN_MASK="${RESULTS_DIR}/centroids_column_mask.npy"

echo "==========================================="
echo "STEP 5: Set up BSCAMPP for phylogenetic placement"
echo "==========================================="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""
echo "Reference: Wedell, Cai & Warnow 2023, IEEE/ACM TCBB 20(2):1417"
echo "  pplacer-BSCAMPP handles backbone trees up to 200,000 leaves."
echo "  Our 167,957-tip tree is used directly — no pruning needed."
echo ""

module purge
module load python3

source "${VENV_DIR}/bin/activate"

# ---------------------------------------------------------------
# Check inputs
# ---------------------------------------------------------------
for f in "$TREE_FILE" "$TREE_LOG" "$ALN_FILE" "$COLUMN_MASK"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        echo "Run steps 3 and 4 first."
        exit 1
    fi
done

TREE_TIPS=$(grep -o ":" "$TREE_FILE" | wc -l)
echo "Reference tree:      $TREE_TIPS tips"
echo "Reference alignment: $ALN_FILE"
echo "FastTree log:        $TREE_LOG"
echo "Column mask:         $COLUMN_MASK"
echo ""

# ---------------------------------------------------------------
# Install BSCAMPP
# ---------------------------------------------------------------
echo "Installing BSCAMPP..."
pip install -q bscampp
python3 -c "import bscampp; print('bscampp OK')"
echo ""

# ---------------------------------------------------------------
# Write per-study placement script
# ---------------------------------------------------------------
PLACEMENT_SCRIPT="${SCRIPT_DIR}/scripts/run_placement.sh"
cat > "$PLACEMENT_SCRIPT" << 'PLACEMENT_EOF'
#!/bin/bash
# Place query 16S sequences into the SILVA 138.2 NR99 bacterial reference tree
# using pplacer-BSCAMPP.
#
# Usage:
#   bash scripts/run_placement.sh query_reads.fasta output_dir output_name [n_cpus]
#
#   $1 = FASTA of unaligned query sequences (ASVs / OTUs / reads)
#   $2 = output directory (will be created if needed)
#   $3 = output name prefix (e.g. "my_study")
#   $4 = number of CPUs (default: 8)
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
VENV_DIR="${SCRIPT_DIR}/venv"

QUERY_FASTA="${1:?Usage: $0 query.fasta output_dir output_name [n_cpus]}"
OUTPUT_DIR="${2:?Usage: $0 query.fasta output_dir output_name [n_cpus]}"
OUTPUT_NAME="${3:?Usage: $0 query.fasta output_dir output_name [n_cpus]}"
N_CPUS="${4:-8}"

TREE_FILE="${RESULTS_DIR}/fasttree/silva_138.2_nr99.nwk"
TREE_LOG="${RESULTS_DIR}/fasttree/silva_138.2_nr99.log"
ALN_FILE="${RESULTS_DIR}/centroids_masked.fasta"
COLUMN_MASK="${RESULTS_DIR}/centroids_column_mask.npy"
QUERY_ALIGNED_RAW="${OUTPUT_DIR}/${OUTPUT_NAME}_aligned_raw.fasta"
QUERY_ALIGNED="${OUTPUT_DIR}/${OUTPUT_NAME}_aligned.fasta"

mkdir -p "$OUTPUT_DIR"
source "${VENV_DIR}/bin/activate"
module load infernal 2>/dev/null || true

# Step 1: align queries to RF00177.cm to get the same 1533 match columns
echo "Aligning queries to RF00177.cm..."
cmalign --cpu "$N_CPUS" \
    --matchonly \
    --dnaout \
    --outformat afa \
    "${SCRIPT_DIR}/models/RF00177.cm" \
    "$QUERY_FASTA" \
    > "$QUERY_ALIGNED_RAW"

# Step 2: apply the same column mask used on the reference alignment.
# Reference and query MUST have identical columns for pplacer.
echo "Applying reference column mask..."
python3 "${SCRIPT_DIR}/scripts/mask_alignment.py" \
    "$QUERY_ALIGNED_RAW" \
    "$QUERY_ALIGNED" \
    --apply-mask "$COLUMN_MASK"

# Step 3: place with pplacer-BSCAMPP
echo "Running pplacer-BSCAMPP..."
run_bscampp.py \
    -i "$TREE_LOG" \
    -t "$TREE_FILE" \
    -a "$ALN_FILE" \
    -q "$QUERY_ALIGNED" \
    -d "$OUTPUT_DIR" \
    -o "$OUTPUT_NAME" \
    --threads "$N_CPUS"

echo ""
echo "Placement complete: ${OUTPUT_DIR}/${OUTPUT_NAME}.jplace"
PLACEMENT_EOF

chmod +x "$PLACEMENT_SCRIPT"
echo "Placement script written: $PLACEMENT_SCRIPT"

deactivate

# ---------------------------------------------------------------
# Log to pipeline.log
# ---------------------------------------------------------------
{
    echo ""
    echo "=== STEP 5: BSCAMPP Setup [$(date)] ==="
    echo "Reference tree:  $TREE_TIPS tips (full tree, no pruning)"
    echo "Model:           FastTree log ($TREE_LOG)"
    echo "Method:          pplacer-BSCAMPP (Wedell et al. 2023, pip install bscampp)"
    echo "Wall time: ${SECONDS} seconds"
} >> "$PIPELINE_LOG"

echo ""
echo "✓✓✓ PIPELINE COMPLETE ✓✓✓"
echo ""
echo "Reference files for placement:"
echo "  Tree:         $TREE_FILE"
echo "  Alignment:    $ALN_FILE"
echo "  FastTree log: $TREE_LOG"
echo "  Column mask:  $COLUMN_MASK"
echo ""
echo "To place queries from a study:"
echo "  bash scripts/run_placement.sh query.fasta results/placements my_study 8"
