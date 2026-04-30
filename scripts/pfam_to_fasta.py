#!/usr/bin/env python3
"""
Convert cmalign Pfam (Stockholm) format alignment to FASTA.

Handles multi-block Pfam files where sequences are split across blocks.

Usage:
  python pfam_to_fasta.py input.pfam output.fasta
"""

import sys
import argparse


def pfam_to_fasta(input_path, output_path):
    sequences = {}
    order = []

    with open(input_path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith('#') or line == '//':
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                name, seq = parts
                if name not in sequences:
                    sequences[name] = []
                    order.append(name)
                sequences[name].append(seq)

    with open(output_path, 'w') as f:
        for name in order:
            f.write(f'>{name}\n{"".join(sequences[name])}\n')

    print(f"Converted {len(sequences)} sequences to {output_path}")
    return len(sequences)


def main():
    parser = argparse.ArgumentParser(
        description='Convert Pfam/Stockholm alignment to FASTA'
    )
    parser.add_argument('input', help='Input Pfam file')
    parser.add_argument('output', help='Output FASTA file')
    args = parser.parse_args()

    count = pfam_to_fasta(args.input, args.output)
    print(f"✓ Done: {count} sequences")
    return 0


if __name__ == '__main__':
    sys.exit(main())
