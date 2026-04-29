#!/bin/bash
#
# Setup script for SILVA tree building on Agate supercomputer
# Run this once via SSH to install all dependencies
#
# Usage:
#   ssh agate.msi.umn.edu
#   cd silva_tree_build
#   bash slurm_setup.sh

set -e  # Exit on error

echo "=========================================="
echo "SILVA Tree Build - Agate Setup"
echo "=========================================="

# Check if we're on Agate
if [[ ! -d /home ]]; then
    echo "ERROR: Does not appear to be an HPC environment"
    exit 1
fi

# Load required modules
echo ""
echo "Loading modules..."
module purge  # Clear any existing modules
module load python3  # Load system Python

# Create a virtual environment for this project
echo ""
echo "Creating Python virtual environment..."
if [ -d "venv" ]; then
    echo "  Virtual environment already exists. Skipping..."
else
    python3 -m venv venv
    echo "  Created: venv/"
fi

# Activate virtual environment
source venv/bin/activate
echo "  Activated: venv/"

# Upgrade pip
echo ""
echo "Upgrading pip..."
pip install --upgrade pip setuptools wheel > /dev/null

# Install Python dependencies
echo ""
echo "Installing Python packages..."
pip install -r requirements.txt

# Check for FastTree (should be available as module)
echo ""
echo "Checking FastTree availability..."
module load fasttree 2>/dev/null && {
    echo "  ✓ FastTree available as module"
    FASTTREE_AVAILABLE=1
} || {
    echo "  ⚠ FastTree not found as module"
    echo "  You may need to:"
    echo "    - Load it manually: module load fasttree"
    echo "    - Or install locally if not available on Agate"
    FASTTREE_AVAILABLE=0
}

echo ""
echo "=========================================="
if [ $FASTTREE_AVAILABLE -eq 1 ]; then
    echo "✓ Setup complete!"
else
    echo "⚠ Setup mostly complete - check FastTree availability"
fi
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Update submit_tree_job.sbatch with:"
echo "     - Path to your SILVA alignment file"
echo "     - Desired CPU/memory allocation"
echo "  2. Submit with: sbatch submit_tree_job.sbatch"
echo ""
