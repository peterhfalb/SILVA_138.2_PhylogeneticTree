#!/bin/bash
#SBATCH --job-name=silva-tree
#SBATCH --partition=agsmall
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=logs/silva-tree-%j.out
#SBATCH --error=logs/silva-tree-%j.err

set -e

echo "======================================"
echo "SILVA 138.2 Tree Build - SLURM Job"
echo "======================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_JOB_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start time: $(date)"
echo "======================================"
echo ""

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
cd "$SCRIPT_DIR"
echo "Working directory: $(pwd)"

mkdir -p logs results

echo "Loading modules..."
module purge
module load python3
module load fasttree

N_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export OMP_NUM_THREADS=$N_THREADS
echo "Threading: OMP_NUM_THREADS=$N_THREADS"

echo "Activating Python environment..."
source "$SCRIPT_DIR/venv/bin/activate"

SILVA_FILE="$SCRIPT_DIR/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz"

echo "Checking input file..."
if [ ! -f "$SILVA_FILE" ]; then
    echo "ERROR: SILVA file not found: $SILVA_FILE"
    exit 1
fi
echo "  ✓ Found: $(basename $SILVA_FILE)"
echo "  Size: $(du -h $SILVA_FILE | cut -f1)"

echo ""
echo "Starting tree building..."
python3 "$SCRIPT_DIR/build_silva_tree_clean.py" \
    "$SILVA_FILE" \
    --output-dir "$SCRIPT_DIR/results" \
    --threads "$N_THREADS"

if [ -f "$SCRIPT_DIR/results/silva-138.2-rooted-tree.nwk" ]; then
    echo ""
    echo "======================================"
    echo "✓ SUCCESS"
    echo "======================================"
    echo "Results saved to: $SCRIPT_DIR/results/"
    echo ""
    echo "Output files:"
    ls -lh "$SCRIPT_DIR/results"/ | grep -E "\.(fasta|nwk|tsv|log)$" | awk '{print "  " $9 " (" $5 ")"}'
    echo ""
    echo "Next steps:"
    echo "  1. Download results:"
    echo "     scp -r agate.msi.umn.edu:/users/9/falb0011/silva_tree_build/SILVA_138/results ~/silva_results"
    echo ""
    echo "  2. Visualize tree:"
    echo "     https://itol.embl.de/ (upload silva-138.2-rooted-tree.nwk)"
else
    echo ""
    echo "======================================"
    echo "✗ ERROR: Tree building failed"
    echo "======================================"
    echo "Check logs: $SCRIPT_DIR/logs/silva-tree-$SLURM_JOB_ID.err"
    exit 1
fi

echo ""
echo "End time: $(date)"
