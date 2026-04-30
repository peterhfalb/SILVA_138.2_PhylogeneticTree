#!/usr/bin/env python3
"""
Mask gappy columns in a FASTA alignment.

Removes columns with gap percentage above threshold, following Ben Kaehler's approach.
Input FASTA can be from cmalign or any pre-aligned sequences.

Usage:
  python mask_alignment.py input.fasta output_masked.fasta [--gap-threshold 99.56]
"""

import sys
import argparse
from pathlib import Path
import numpy as np
import skbio


def load_alignment(filepath):
    """Load aligned sequences from FASTA."""
    print(f"Loading alignment from {filepath}...")
    sequences = []

    for i, seq in enumerate(skbio.io.read(filepath, format='fasta', constructor=skbio.DNA)):
        sequences.append(seq)
        if (i + 1) % 50000 == 0:
            print(f"  Loaded {i + 1} sequences")

    print(f"Total sequences: {len(sequences)}")
    return sequences


def mask_gappy_columns(sequences, gap_threshold_pct=99.56):
    """
    Identify and filter out gappy columns.

    Keeps columns with gap% < threshold.

    Args:
        sequences: List of aligned DNA sequences
        gap_threshold_pct: Gap percentage threshold (default 99.56%)

    Returns:
        Boolean mask array for columns to keep
    """
    print("\nAnalyzing gap patterns...")

    seq_array = np.array([seq.values for seq in sequences])
    n_seqs, n_cols = seq_array.shape
    print(f"Alignment dimensions: {n_seqs} sequences × {n_cols} columns")

    num_gaps = np.zeros(n_cols, dtype=int)
    for j in range(n_cols):
        num_gaps[j] = (seq_array[:, j] == b'.').sum() + (seq_array[:, j] == b'-').sum()
        if (j + 1) % 10000 == 0:
            print(f"  Analyzed {j + 1} columns")

    num_ok = int((gap_threshold_pct / 100.0) * n_seqs)
    keep_mask = num_gaps <= (n_seqs - num_ok)

    n_kept = keep_mask.sum()

    print(f"\nGap masking results:")
    print(f"  Original: {n_cols} columns")
    print(f"  Retained: {n_kept} columns")
    print(f"  Removed: {n_cols - n_kept} gappy columns")
    print(f"  Threshold: {gap_threshold_pct:.2f}%")

    return keep_mask


def export_masked_alignment(sequences, keep_mask, output_path):
    """Export masked alignment to FASTA."""
    print(f"\nExporting masked alignment to {output_path}...")

    seq_count = 0
    with open(output_path, 'w') as f:
        for seq in sequences:
            masked_values = seq.values[keep_mask]
            seq_str = ''.join([b.decode('utf-8') if isinstance(b, bytes) else b
                              for b in masked_values])
            seq_id = seq.metadata.get('id', f'seq_{seq_count}')
            f.write(f">{seq_id}\n{seq_str}\n")
            seq_count += 1

            if (seq_count % 50000) == 0:
                print(f"  Exported {seq_count} sequences")

    print(f"Total exported: {seq_count} sequences")


def main():
    parser = argparse.ArgumentParser(
        description='Mask gappy columns in aligned FASTA'
    )
    parser.add_argument('input', help='Input aligned FASTA')
    parser.add_argument('output', help='Output masked FASTA')
    parser.add_argument('--gap-threshold', type=float, default=99.56,
                       help='Gap percentage threshold (default 99.56)')

    args = parser.parse_args()

    sequences = load_alignment(args.input)
    keep_mask = mask_gappy_columns(sequences, args.gap_threshold)
    export_masked_alignment(sequences, keep_mask, args.output)

    print("\n✓ Masking complete")
    return 0


if __name__ == '__main__':
    sys.exit(main())
