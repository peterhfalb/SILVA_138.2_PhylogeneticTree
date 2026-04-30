#!/usr/bin/env python3
"""
Filter centroid sequences to bacteria only, using original SILVA taxonomy.

VSEARCH strips descriptions from headers, so we look up taxonomy from the
original SILVA FASTA by matching accession numbers.

Usage:
  python filter_bacteria.py centroids.fasta silva.fasta.gz output_bacteria.fasta
"""

import sys
import gzip
import argparse


def load_centroid_ids(centroid_fasta):
    """Load set of centroid accession IDs."""
    ids = set()
    with open(centroid_fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line.strip()[1:].split()[0])
    print(f"Centroid IDs loaded: {len(ids)}")
    return ids


def build_taxonomy_map(silva_file, centroid_ids):
    """
    Scan original SILVA FASTA to get taxonomy for centroid IDs.
    Returns dict of {seqid: full_taxonomy_string}.
    """
    taxonomy = {}
    opener = gzip.open if silva_file.endswith(".gz") else open

    with opener(silva_file, "rt") as f:
        for line in f:
            if line.startswith(">"):
                parts = line.strip()[1:].split(" ", 1)
                seq_id = parts[0]
                if seq_id in centroid_ids:
                    taxonomy[seq_id] = parts[1] if len(parts) > 1 else ""

    print(f"Taxonomy found for: {len(taxonomy)} / {len(centroid_ids)} centroids")
    return taxonomy


def filter_bacteria(centroid_fasta, taxonomy_map, output_fasta):
    """Write bacteria-only centroids to output FASTA."""
    kept = 0
    skipped = 0
    no_taxonomy = 0

    with open(centroid_fasta) as fin, open(output_fasta, "w") as fout:
        write_seq = False
        current_id = None
        for line in fin:
            if line.startswith(">"):
                seq_id = line.strip()[1:].split()[0]
                tax = taxonomy_map.get(seq_id, "")
                if not tax:
                    no_taxonomy += 1
                    write_seq = False
                elif tax.startswith("Bacteria;"):
                    write_seq = True
                    # Write header with taxonomy restored
                    fout.write(f">{seq_id} {tax}\n")
                    kept += 1
                else:
                    write_seq = False
                    skipped += 1
            elif write_seq:
                fout.write(line)

    print(f"Bacteria kept:          {kept}")
    print(f"Non-bacteria skipped:   {skipped}")
    print(f"No taxonomy found:      {no_taxonomy}")
    return kept


def main():
    parser = argparse.ArgumentParser(
        description="Filter centroids to bacteria only using original SILVA taxonomy"
    )
    parser.add_argument("centroids", help="Centroid FASTA from VSEARCH")
    parser.add_argument("silva", help="Original SILVA FASTA (gzipped or plain)")
    parser.add_argument("output", help="Output FASTA (bacteria only)")
    args = parser.parse_args()

    centroid_ids = load_centroid_ids(args.centroids)
    taxonomy_map = build_taxonomy_map(args.silva, centroid_ids)
    kept = filter_bacteria(args.centroids, taxonomy_map, args.output)

    print(f"\n✓ Done. {kept} bacteria centroids written to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
