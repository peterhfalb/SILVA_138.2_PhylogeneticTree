#!/usr/bin/env python3
"""
Build a phylogenetic tree from SILVA 138.2 NR99 using FastTree.

This implements Ben Kaehler's approach from:
https://github.com/BenKaehler/silva-trees

Steps:
  1. Load SILVA 138.2 SSURef NR99 aligned sequences (RNA)
  2. Convert RNA to DNA
  3. Mask gappy columns (keep <99.56% gaps)
  4. Export masked alignment to FASTA
  5. Run FastTree and capture detailed logs
  6. Extract taxonomy from sequence headers
  7. Root tree at midpoint

Outputs:
  - silva-138.2-masked-aln.fasta: Masked alignment
  - silva-138.2-tree.nwk: Unrooted tree (Newick format)
  - silva-138.2-rooted-tree.nwk: Midpoint-rooted tree
  - silva-138.2-fasttree.log: FastTree output with model parameters
  - silva-138.2-taxonomy.tsv: Taxonomy table (ID and taxonomy)
"""

import sys
import gzip
import subprocess
from pathlib import Path
import numpy as np
import skbio


def load_and_convert_alignment(filepath, count_every=50000):
    """
    Load SILVA RNA alignment and convert to DNA.

    Args:
        filepath: Path to SILVA FASTA file (gzipped)
        count_every: Print progress every N sequences

    Returns:
        list: Converted DNA sequences
    """
    print("Loading SILVA alignment (RNA -> DNA conversion)...")
    sequences = []

    # Handle gzipped files
    if filepath.endswith('.gz'):
        handle = gzip.open(filepath, 'rt')
    else:
        handle = open(filepath, 'r')

    with handle:
        for i, rna_seq in enumerate(skbio.io.read(handle, format='fasta', constructor=skbio.RNA)):
            # Convert RNA to DNA (U -> T)
            dna_seq = rna_seq.reverse_transcribe()
            sequences.append(dna_seq)

            if (i + 1) % count_every == 0:
                print(f"  Loaded {i + 1} sequences")

    print(f"Total sequences loaded: {len(sequences)}")
    return sequences


def mask_gappy_columns(sequences, gap_threshold_pct=99.56, report_every=10000):
    """
    Identify and filter out gappy columns.

    Follows Ben Kaehler's approach: keep columns with gap% < threshold.

    Args:
        sequences: List of aligned DNA sequences
        gap_threshold_pct: Gap percentage threshold (default 99.56%)
        report_every: Print progress every N columns

    Returns:
        tuple: (boolean mask array, filtered sequences)
    """
    print("\nAnalyzing gap patterns...")

    # Convert to numpy array for column access
    seq_array = np.array([seq.values for seq in sequences])
    n_seqs, n_cols = seq_array.shape

    # Count gaps in each column
    num_gaps = np.zeros(n_cols, dtype=int)
    for j in range(n_cols):
        num_gaps[j] = (seq_array[:, j] == b'.').sum() + (seq_array[:, j] == b'-').sum()
        if (j + 1) % report_every == 0:
            print(f"  Analyzed {j + 1} columns")

    # Determine which columns to keep
    # Using the approach from Ben's notebook:
    # num_ok = 1907 (number of sequences with gaps)
    # Keep columns where gaps <= n_seqs - num_ok
    num_ok = int((gap_threshold_pct / 100.0) * n_seqs)
    keep_mask = num_gaps <= (n_seqs - num_ok)

    n_kept = keep_mask.sum()
    actual_gap_pct = 100.0 * (num_gaps[keep_mask].max()) / n_seqs

    print(f"\nGap masking results:")
    print(f"  Original alignment: {n_cols} columns")
    print(f"  Keeping: {n_kept} columns")
    print(f"  Removing: {n_cols - n_kept} gappy columns")
    print(f"  Gap threshold: {gap_threshold_pct:.2f}%")

    return keep_mask, sequences


def export_masked_alignment(sequences, keep_mask, output_path):
    """
    Export masked alignment to proper FASTA format.

    Args:
        sequences: List of DNA sequences
        keep_mask: Boolean array of columns to keep
        output_path: Path to output FASTA file
    """
    print(f"\nExporting masked alignment to {output_path}...")

    seq_count = 0
    with open(output_path, 'w') as f:
        for seq in sequences:
            # Extract masked columns
            masked_values = seq.values[keep_mask]

            # Convert bytes to string
            if isinstance(masked_values, np.ndarray):
                seq_str = ''.join([b.decode('utf-8') if isinstance(b, bytes) else b
                                  for b in masked_values])
            else:
                seq_str = masked_values.decode('utf-8') if isinstance(masked_values, bytes) else str(masked_values)

            # Get sequence ID
            seq_id = seq.metadata.get('id', f'seq_{seq_count}')

            # Write FASTA format
            f.write(f">{seq_id}\n{seq_str}\n")
            seq_count += 1

            if (seq_count % 50000) == 0:
                print(f"  Exported {seq_count} sequences")

    print(f"  Total: {seq_count} sequences exported")


def extract_taxonomy(sequences, output_path):
    """
    Extract taxonomy information from SILVA sequence headers.

    SILVA header format: ID description;taxonomy

    Args:
        sequences: List of sequences with metadata
        output_path: Path to output TSV file
    """
    print(f"\nExtracting taxonomy to {output_path}...")

    with open(output_path, 'w') as f:
        f.write("seqid\ttaxonomy\n")

        for seq in sequences:
            seq_id = seq.metadata.get('id', '')
            description = seq.metadata.get('description', '')

            # Try to extract taxonomy from description
            # SILVA format: full header with taxonomy at the end after semicolon
            if description and ';' in description:
                # Get everything after the last semicolon as taxonomy
                taxonomy = description.split(';')[-1].strip()
            else:
                # Fallback: use the description if available
                taxonomy = description.strip() if description else 'unknown'

            f.write(f"{seq_id}\t{taxonomy}\n")

    print(f"  Taxonomy table saved")


def run_fasttree(alignment_fasta, output_tree, output_log, n_threads=4):
    """
    Run FastTree to build phylogenetic tree.

    Args:
        alignment_fasta: Path to FASTA alignment
        output_tree: Path to output Newick tree
        output_log: Path to save FastTree log/output
        n_threads: Number of threads (via OMP_NUM_THREADS)

    Returns:
        bool: True if successful
    """
    print(f"\nBuilding tree with FastTree ({n_threads} threads)...")

    # Try to find FastTree binary
    fasttree_options = ['FastTreeMP', 'FastTree']
    fasttree_cmd = None

    for cmd in fasttree_options:
        try:
            subprocess.run([cmd, '-help'], capture_output=True, timeout=5, check=False)
            fasttree_cmd = cmd
            print(f"  Using: {fasttree_cmd}")
            break
        except FileNotFoundError:
            continue

    if not fasttree_cmd:
        print(f"ERROR: FastTree not found")
        return False

    # Build command
    cmd = [
        fasttree_cmd,
        '-quote',        # Handle special chars in names
        '-nt',           # Nucleotide sequences
        alignment_fasta
    ]

    # Run FastTree and capture output
    try:
        print(f"  Running: {' '.join(cmd)}")
        with open(output_log, 'w') as log_file:
            with open(output_tree, 'w') as tree_file:
                result = subprocess.run(
                    cmd,
                    stdout=tree_file,
                    stderr=log_file,
                    text=True
                )

        if result.returncode != 0:
            print(f"ERROR: FastTree exited with code {result.returncode}")
            with open(output_log, 'r') as f:
                print(f.read())
            return False

        print(f"  Tree saved to {output_tree}")
        print(f"  Logs saved to {output_log}")

        # Print FastTree parameters from log
        print(f"\nFastTree parameters:")
        with open(output_log, 'r') as f:
            for line in f:
                if any(key in line for key in ['FastTree', 'threads', 'Model', 'Nucleotide', 'Support']):
                    print(f"  {line.rstrip()}")

        return True

    except Exception as e:
        print(f"ERROR: {e}")
        return False


def midpoint_root(unrooted_tree, rooted_tree):
    """
    Root tree at its midpoint.

    Args:
        unrooted_tree: Path to unrooted Newick tree
        rooted_tree: Path to save rooted tree

    Returns:
        bool: True if successful
    """
    print(f"\nMidpoint rooting tree...")

    try:
        with open(unrooted_tree, 'r') as f:
            tree_str = f.read().strip()

        tree = skbio.tree.TreeNode.read(tree_str, format='newick')
        tree_rooted = tree.root_at_midpoint()

        with open(rooted_tree, 'w') as f:
            f.write(str(tree_rooted) + ';')

        print(f"  Rooted tree saved to {rooted_tree}")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        return False


def main():
    """Main workflow."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Build SILVA 138.2 phylogenetic tree (Ben Kaehler approach)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python build_silva_tree_clean.py SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz
  python build_silva_tree_clean.py alignment.fasta.gz -o output/
  python build_silva_tree_clean.py alignment.fasta.gz -t 16
        """
    )

    parser.add_argument('alignment', help='SILVA FASTA alignment (gzipped)')
    parser.add_argument('-o', '--output-dir', default='.', help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define output files
    masked_fasta = output_dir / 'silva-138.2-masked-aln.fasta'
    unrooted_tree = output_dir / 'silva-138.2-tree.nwk'
    rooted_tree = output_dir / 'silva-138.2-rooted-tree.nwk'
    tree_log = output_dir / 'silva-138.2-fasttree.log'
    taxonomy_file = output_dir / 'silva-138.2-taxonomy.tsv'

    # === STEP 1: Load and convert alignment ===
    sequences = load_and_convert_alignment(args.alignment)
    if not sequences:
        print("ERROR: No sequences loaded")
        return 1

    # === STEP 2: Mask gappy columns ===
    keep_mask, sequences = mask_gappy_columns(sequences)

    # === STEP 3: Export masked alignment ===
    export_masked_alignment(sequences, keep_mask, str(masked_fasta))

    # === STEP 4: Extract taxonomy ===
    extract_taxonomy(sequences, str(taxonomy_file))

    # === STEP 5: Build tree with FastTree ===
    if not run_fasttree(str(masked_fasta), str(unrooted_tree), str(tree_log), n_threads=args.threads):
        return 1

    # === STEP 6: Root tree at midpoint ===
    if not midpoint_root(str(unrooted_tree), str(rooted_tree)):
        return 1

    # === Summary ===
    print("\n" + "=" * 60)
    print("✓ SUCCESS: Tree building complete!")
    print("=" * 60)
    print("\nOutput files:")
    print(f"  Alignment: {masked_fasta}")
    print(f"  Unrooted tree: {unrooted_tree}")
    print(f"  Rooted tree: {rooted_tree}")
    print(f"  FastTree log: {tree_log}")
    print(f"  Taxonomy: {taxonomy_file}")
    print("\nVisualize with:")
    print(f"  iTOL: https://itol.embl.de/ (upload {rooted_tree})")

    return 0


if __name__ == '__main__':
    sys.exit(main())
