#!/usr/bin/env python3
"""
Extract full SILVA taxonomy lineage from sequence headers.

SILVA header format: >ID full_semicolon_delimited_lineage

Outputs a TSV with seqid and full lineage for all sequences.

Usage:
  python extract_taxonomy.py centroids.fasta output_taxonomy.tsv
"""

import sys
import argparse


def extract_taxonomy(fasta_path, output_path):
    """
    Extract full SILVA lineage from sequence headers.

    SILVA format: >ID taxonomy_string (entire semicolon-delimited lineage)
    Parses headers directly to avoid skbio DNA/RNA character validation.

    Args:
        fasta_path: Path to input FASTA
        output_path: Path to output TSV
    """
    print(f"Reading {fasta_path}...")

    seq_count = 0

    with open(fasta_path) as in_f, open(output_path, 'w') as out_f:
        out_f.write("seqid\ttaxonomy\n")

        for line in in_f:
            if not line.startswith('>'):
                continue
            parts = line[1:].rstrip().split(' ', 1)
            seq_id = parts[0]
            taxonomy = parts[1].strip() if len(parts) > 1 else 'unknown'
            out_f.write(f"{seq_id}\t{taxonomy}\n")
            seq_count += 1

            if seq_count % 50000 == 0:
                print(f"  Processed {seq_count} sequences")

    print(f"✓ Extracted {seq_count} sequences to {output_path}")
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='Extract full SILVA taxonomy lineage from sequence headers'
    )
    parser.add_argument('fasta', help='Input FASTA file')
    parser.add_argument('output', help='Output TSV file (seqid + full lineage)')

    args = parser.parse_args()
    return extract_taxonomy(args.fasta, args.output)


if __name__ == '__main__':
    sys.exit(main())
