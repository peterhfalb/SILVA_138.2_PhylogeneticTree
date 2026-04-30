#!/bin/bash
# Download SILVA NR99 138.2 unaligned sequences
# Run interactively (not via SLURM): bash scripts/01_download_silva.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SILVA_FILE="SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
SILVA_URL="https://ftp.arb-silva.de/release_138.2/Exports/${SILVA_FILE}"
LOGS_DIR="${SCRIPT_DIR}/logs"

mkdir -p "$LOGS_DIR"

echo "=== Downloading SILVA NR99 138.2 unaligned sequences ==="
echo "Source: $SILVA_URL"
echo "Destination: ${SCRIPT_DIR}/${SILVA_FILE}"
echo ""

if [ -f "${SCRIPT_DIR}/${SILVA_FILE}" ]; then
    echo "✓ File already exists: ${SILVA_FILE}"
    ls -lh "${SCRIPT_DIR}/${SILVA_FILE}"
else
    echo "Downloading (this may take a few minutes)..."
    cd "$SCRIPT_DIR"
    wget -c "$SILVA_URL"
    echo ""
    echo "✓ Download complete"
    ls -lh "${SILVA_FILE}"
fi

echo ""
echo "Next step: Submit clustering job"
echo "  sbatch scripts/02_cluster.sbatch"
