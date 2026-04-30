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
import skbio


def extract_taxonomy(fasta_path, output_path):
    """
    Extract full SILVA lineage from sequence headers.

    SILVA format: >ID taxonomy_string (entire semicolon-delimited lineage)
    This extracts the FULL lineage, not just the last element.

    Args:
        fasta_path: Path to input FASTA
        output_path: Path to output TSV
    """
    print(f"Reading {fasta_path}...")

    seq_count = 0
    written_count = 0

    with open(output_path, 'w') as out_f:
        # Write header
        out_f.write("seqid\ttaxonomy\n")

        # Read sequences
        for seq in skbio.io.read(fasta_path, format='fasta', constructor=skbio.DNA):
            seq_id = seq.metadata.get('id', '')
            description = seq.metadata.get('description', '')

            # In scikit-bio, 'description' is everything after the first space in the header
            # For SILVA format: >ID TAXONOMY
            # So description = TAXONOMY (the full lineage string)
            if description:
                taxonomy = description.strip()
            else:
                taxonomy = 'unknown'

            out_f.write(f"{seq_id}\t{taxonomy}\n")
            written_count += 1
            seq_count += 1

            if seq_count % 50000 == 0:
                print(f"  Processed {seq_count} sequences")

    print(f"✓ Extracted {written_count} sequences to {output_path}")
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
