#!/bin/bash
# Download SILVA NR99 138.2 unaligned sequences and taxonomy files.
# Run interactively on the login node (not via SLURM): bash scripts/01_download_silva.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOGS_DIR="${SCRIPT_DIR}/logs"

SILVA_FILE="SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
SILVA_URL="https://ftp.arb-silva.de/release_138.2/Exports/${SILVA_FILE}"

TAXMAP_GZ="taxmap_slv_ssu_ref_nr_138.2.txt.gz"
TAXMAP_FILE="taxmap_slv_ssu_ref_nr_138.2.txt"
TAXMAP_URL="https://ftp.arb-silva.de/release_138.2/Exports/taxonomy/${TAXMAP_GZ}"

mkdir -p "$LOGS_DIR"

# ---------------------------------------------------------------
# SILVA sequences
# ---------------------------------------------------------------
echo "=== Downloading SILVA NR99 138.2 unaligned sequences ==="
echo "Source: $SILVA_URL"
echo ""

if [ -f "${SCRIPT_DIR}/${SILVA_FILE}" ]; then
    echo "✓ Already exists: ${SILVA_FILE}"
    ls -lh "${SCRIPT_DIR}/${SILVA_FILE}"
else
    echo "Downloading (this may take a few minutes)..."
    cd "$SCRIPT_DIR"
    wget -c "$SILVA_URL"
    echo "✓ Download complete"
    ls -lh "${SILVA_FILE}"
fi

# ---------------------------------------------------------------
# SILVA taxmap (provides NCBI taxonomy IDs per sequence)
# ---------------------------------------------------------------
echo ""
echo "=== Downloading SILVA taxmap (NCBI taxonomy IDs) ==="
echo "Source: $TAXMAP_URL"
echo ""

if [ -f "${SCRIPT_DIR}/${TAXMAP_FILE}" ]; then
    echo "✓ Already exists: ${TAXMAP_FILE}"
    ls -lh "${SCRIPT_DIR}/${TAXMAP_FILE}"
elif [ -f "${SCRIPT_DIR}/${TAXMAP_GZ}" ]; then
    echo "Decompressing existing ${TAXMAP_GZ}..."
    gunzip -k "${SCRIPT_DIR}/${TAXMAP_GZ}"
    echo "✓ Done"
    ls -lh "${SCRIPT_DIR}/${TAXMAP_FILE}"
else
    cd "$SCRIPT_DIR"
    wget -c "$TAXMAP_URL"
    echo "Decompressing..."
    gunzip -k "${TAXMAP_GZ}"
    echo "✓ Download and decompress complete"
    ls -lh "${TAXMAP_FILE}"
fi

# ---------------------------------------------------------------
# NCBI taxonomy database (required for taxit taxtable in step 5)
# This must be created here (login node) because compute nodes may
# not have outbound internet access.
# ---------------------------------------------------------------
echo ""
echo "=== Building NCBI taxonomy database ==="
TAXONOMY_DIR="${SCRIPT_DIR}/results/taxonomy"
TAXONOMY_DB="${TAXONOMY_DIR}/taxonomy.db"

mkdir -p "$TAXONOMY_DIR"

if [ -f "$TAXONOMY_DB" ]; then
    echo "✓ taxonomy.db already exists: $TAXONOMY_DB"
    ls -lh "$TAXONOMY_DB"
else
    if ! command -v taxit &>/dev/null; then
        echo "NOTE: taxit not found on PATH."
        echo "  Set up the venv first, then re-run this script:"
        echo "    module load python3"
        echo "    python3 -m venv venv"
        echo "    source venv/bin/activate"
        echo "    pip install -r requirements.txt"
        echo "  Then: source venv/bin/activate && bash scripts/01_download_silva.sh"
        echo ""
        echo "Skipping taxonomy.db creation — run 'taxit new_database $TAXONOMY_DB' manually."
    else
        echo "Running taxit new_database (downloads NCBI taxdmp.zip, may take 5–10 min)..."
        cd "$TAXONOMY_DIR"
        taxit new_database "$TAXONOMY_DB"
        cd "$SCRIPT_DIR"
        echo "✓ taxonomy.db created"
        ls -lh "$TAXONOMY_DB"
    fi
fi

echo ""
echo "Next step: Submit clustering job"
echo "  sbatch scripts/02_cluster.sbatch"
