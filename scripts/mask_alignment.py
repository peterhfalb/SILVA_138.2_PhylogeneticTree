#!/usr/bin/env python3
"""
Mask gappy columns in a FASTA alignment.

Removes columns with gap percentage above threshold, following Ben Kaehler's approach.
Input FASTA can be from cmalign or any pre-aligned sequences.

Two modes:
  1. Derive mask from reference alignment (use during tree building):
       python mask_alignment.py input.fasta output.fasta [--gap-threshold 99.56] [--mask-out mask.npy]
     --mask-out saves the retained column indices so the same columns can be
     applied to query alignments later.

  2. Apply a saved mask to a query alignment (use during placement):
       python mask_alignment.py query.fasta output.fasta --apply-mask mask.npy
     The query alignment must have the same total number of columns as the
     original reference alignment before masking (i.e. 1533 cmalign match columns).
"""

import sys
import argparse
from pathlib import Path
import numpy as np
import skbio


def load_alignment(filepath):
    """Load aligned sequences from FASTA, handling both DNA and RNA input."""
    print(f"Loading alignment from {filepath}...")
    sequences = []

    # Try DNA first; fall back to RNA and convert (cmalign outputs RNA by default)
    try:
        for i, seq in enumerate(skbio.io.read(filepath, format='fasta', constructor=skbio.DNA)):
            sequences.append(seq)
            if (i + 1) % 50000 == 0:
                print(f"  Loaded {i + 1} sequences")
    except ValueError:
        print("  RNA characters detected, converting U→T...")
        sequences = []
        for i, seq in enumerate(skbio.io.read(filepath, format='fasta', constructor=skbio.RNA)):
            sequences.append(seq.reverse_transcribe())
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


def load_mask(mask_path):
    """Load a previously saved column mask (boolean array) from a .npy file."""
    mask = np.load(mask_path)
    print(f"Loaded column mask: {mask.sum()} of {len(mask)} columns retained ({mask_path})")
    return mask


def main():
    parser = argparse.ArgumentParser(
        description='Mask gappy columns in aligned FASTA'
    )
    parser.add_argument('input', help='Input aligned FASTA')
    parser.add_argument('output', help='Output masked FASTA')
    parser.add_argument('--gap-threshold', type=float, default=99.56,
                       help='Gap percentage threshold for deriving mask (default 99.56)')
    parser.add_argument('--mask-out', metavar='FILE',
                       help='Save retained column indices to this .npy file '
                            '(use when masking the reference alignment)')
    parser.add_argument('--apply-mask', metavar='FILE',
                       help='Apply a previously saved column mask instead of '
                            'deriving one from gap frequency '
                            '(use when masking query alignments for placement)')

    args = parser.parse_args()

    sequences = load_alignment(args.input)

    if args.apply_mask:
        keep_mask = load_mask(args.apply_mask)
        n_cols = len(sequences[0].values)
        if len(keep_mask) != n_cols:
            print(f"ERROR: mask has {len(keep_mask)} columns but alignment has {n_cols}")
            print("Query alignment must have the same number of columns as the "
                  "original reference alignment before masking.")
            return 1
    else:
        keep_mask = mask_gappy_columns(sequences, args.gap_threshold)
        if args.mask_out:
            np.save(args.mask_out, keep_mask)
            print(f"Column mask saved to {args.mask_out}")

    export_masked_alignment(sequences, keep_mask, args.output)

    print("\n✓ Masking complete")
    return 0


if __name__ == '__main__':
    sys.exit(main())
