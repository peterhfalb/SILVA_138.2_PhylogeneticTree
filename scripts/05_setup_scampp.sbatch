#!/bin/bash
#SBATCH --job-name=silva-scampp-setup
#SBATCH --partition=agsmall
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=2:00:00
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

mkdir -p "${RESULTS_DIR}/scampp" "$LOGS_DIR"

PIPELINE_LOG="${LOGS_DIR}/pipeline.log"
TREE_FILE="${RESULTS_DIR}/fasttree/silva_138.2_nr99.nwk"
TREE_LOG="${RESULTS_DIR}/fasttree/silva_138.2_nr99.log"
ALN_FILE="${RESULTS_DIR}/centroids_masked.fasta"
REFPKG="${RESULTS_DIR}/scampp/silva_138.2_nr99.refpkg"

echo "==========================================="
echo "STEP 5: Set up SCAMPP for phylogenetic placement"
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
if [ ! -f "$TREE_FILE" ]; then
    echo "ERROR: Tree file not found: $TREE_FILE"
    echo "Run step 4 first: sbatch scripts/04_build_tree.sbatch"
    exit 1
fi

if [ ! -f "$ALN_FILE" ]; then
    echo "ERROR: Masked alignment not found: $ALN_FILE"
    exit 1
fi

TREE_TIPS=$(grep -o ":" "$TREE_FILE" | wc -l)
echo "Backbone tree: $TREE_TIPS tips"
echo ""

# ---------------------------------------------------------------
# Install BSCAMPP and dependencies
# ---------------------------------------------------------------
echo "Installing BSCAMPP and dependencies..."
pip install bscampp taxtastic
echo ""

# Verify bscampp is importable
python3 -c "import bscampp; print('bscampp OK')"

# ---------------------------------------------------------------
# Locate or download pplacer binary
# ---------------------------------------------------------------
#echo "Locating pplacer..."
#PPLACER_BIN=""

# Try module, then PATH, then download
#if module spider pplacer 2>/dev/null | grep -qi pplacer; then
#    module load pplacer
#    PPLACER_BIN=$(which pplacer)
#    echo "  Found via module: $PPLACER_BIN"
#elif command -v pplacer &>/dev/null; then
#    PPLACER_BIN=$(which pplacer)
#    echo "  Found in PATH: $PPLACER_BIN"
#else
#   PPLACER_BIN="${SCRIPT_DIR}/bin/pplacer"
#    mkdir -p "${SCRIPT_DIR}/bin"
#   if [ ! -f "$PPLACER_BIN" ]; then
#        echo "  Downloading pplacer v1.1.alpha19 (Linux x86_64)..."
#        wget -q -O /tmp/pplacer.zip \
#            "https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip"
#        unzip -q /tmp/pplacer.zip -d /tmp/pplacer_tmp
#        cp /tmp/pplacer_tmp/*/pplacer "${SCRIPT_DIR}/bin/"
#        cp /tmp/pplacer_tmp/*/guppy   "${SCRIPT_DIR}/bin/" 2>/dev/null || true
#        cp /tmp/pplacer_tmp/*/rppr    "${SCRIPT_DIR}/bin/" 2>/dev/null || true
#        chmod +x "${SCRIPT_DIR}/bin/"*
#        rm -rf /tmp/pplacer.zip /tmp/pplacer_tmp
#    fi
#    echo "  Using local binary: $PPLACER_BIN"
#fi
#
#if [ ! -x "$PPLACER_BIN" ]; then
#    echo "ERROR: pplacer not found or not executable: $PPLACER_BIN"
#    echo "  Download from: https://github.com/matsen/pplacer/releases"
#    exit 1
#fi
#echo "  pplacer: $(${PPLACER_BIN} --version 2>&1 | head -1)"
#echo ""

# ---------------------------------------------------------------
# Create Taxtastic reference package
# FastTree log provides GTR rates + gamma parameters for pplacer.
# Per SCAMPP paper: FastTree parameters work directly with pplacer-SCAMPP.
# ---------------------------------------------------------------
echo "Creating Taxtastic reference package..."
echo "  Tree:   $TREE_FILE"
echo "  Log:    $TREE_LOG"
echo "  Align:  $ALN_FILE"
echo "  Output: $REFPKG"

[ -d "$REFPKG" ] && rm -rf "$REFPKG"

taxit create \
    -l "silva_138.2_nr99_bacteria" \
    -P "$REFPKG" \
    --aln-fasta "$ALN_FILE" \
    --tree "$TREE_FILE" \
    --tree-stats "$TREE_LOG"

echo ""
echo "Reference package contents:"
ls -lh "$REFPKG"/

# ---------------------------------------------------------------
# Write per-study placement script template
# ---------------------------------------------------------------
PLACEMENT_SCRIPT="${SCRIPT_DIR}/scripts/run_placement.sh"
cat > "$PLACEMENT_SCRIPT" << 'PLACEMENT_EOF'
#!/bin/bash
# Place query 16S reads into the SILVA reference tree with pplacer-BSCAMPP.
#
# Usage:
#   bash scripts/run_placement.sh query_reads.fasta output_prefix [n_cpus]
#
#   $1 = FASTA of unaligned query sequences (ASVs / OTUs / reads)
#   $2 = output prefix (e.g. "results/scampp/my_study")
#   $3 = number of CPUs (default: 8)
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
VENV_DIR="${SCRIPT_DIR}/venv"

QUERY_FASTA="${1:?Usage: $0 query.fasta output_prefix [n_cpus]}"
OUTPUT_PREFIX="${2:?Usage: $0 query.fasta output_prefix [n_cpus]}"
N_CPUS="${3:-8}"

ALN_FILE="${RESULTS_DIR}/centroids_masked.fasta"
TREE_FILE="${RESULTS_DIR}/fasttree/silva_138.2_nr99.nwk"
REFPKG="${RESULTS_DIR}/scampp/silva_138.2_nr99.refpkg"
QUERY_ALIGNED="${OUTPUT_PREFIX}_aligned.fasta"
JPLACE_OUT="${OUTPUT_PREFIX}.jplace"

mkdir -p "$(dirname "$OUTPUT_PREFIX")"
source "${VENV_DIR}/bin/activate"
module load infernal 2>/dev/null || true

# Step 1: align queries to RF00177.cm (same model as reference alignment)
echo "Aligning queries to RF00177.cm..."
cmalign --cpu "$N_CPUS" \
    --matchonly \
    --dnaout \
    --outformat afa \
    "${SCRIPT_DIR}/models/RF00177.cm" \
    "$QUERY_FASTA" \
    > "$QUERY_ALIGNED"

# Step 2: place with pplacer-BSCAMPP (B=2000, per Wedell et al. 2023)
echo "Running pplacer-BSCAMPP (B=2000)..."
run_bscampp.py \
    -t "$TREE_FILE" \
    -a "$ALN_FILE" \
    -q "$QUERY_ALIGNED" \
    -d "$REFPKG" \
    -o "$JPLACE_OUT" \
    -b 2000 \
    -j "$N_CPUS"

echo ""
echo "Placement complete: $JPLACE_OUT"
echo "Summarize with guppy: guppy fat $JPLACE_OUT"
PLACEMENT_EOF

chmod +x "$PLACEMENT_SCRIPT"

deactivate

# ---------------------------------------------------------------
# Log to pipeline.log
# ---------------------------------------------------------------
{
    echo ""
    echo "=== STEP 5: SCAMPP Setup [$(date)] ==="
    echo "Backbone tree: $TREE_TIPS tips (full tree, no pruning)"
    echo "Method: pplacer-BSCAMPP (Wedell et al. 2023, pip install bscampp)"
    echo "  B=2000 (placement subtree size, per paper recommendation)"
    echo "Reference package: $REFPKG"
    echo "pplacer binary: $PPLACER_BIN"
    echo "Wall time: ${SECONDS} seconds"
} >> "$PIPELINE_LOG"

echo ""
echo "✓✓✓ PIPELINE COMPLETE ✓✓✓"
echo ""
echo "Summary of output files:"
echo "  - Masked alignment:  ${ALN_FILE}"
echo "  - Full tree:         ${TREE_FILE}"
echo "  - FastTree log:      ${TREE_LOG}"
echo "  - Taxonomy TSV:      ${RESULTS_DIR}/centroids_taxonomy.tsv"
echo "  - SCAMPP refpkg:     ${REFPKG}"
echo ""
echo "To place queries from a study:"
echo "  bash scripts/run_placement.sh /path/to/query_reads.fasta results/scampp/my_study 8"
