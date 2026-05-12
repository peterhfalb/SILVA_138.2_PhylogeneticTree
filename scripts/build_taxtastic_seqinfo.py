#!/usr/bin/env python3
"""
Build seq_info.csv and tax_ids.txt for Taxtastic using the SILVA taxmap file.

The SILVA taxmap provides NCBI taxonomy IDs for each sequence. seq_info.csv and
tax_ids.txt are then used with 'taxit taxtable' to extract the minimal NCBI
taxonomy subtree for the reference package.

SILVA taxmap format (tab-delimited, with header):
  primaryAccession  start  stop  path  organism_name  ncbi_taxid

The seqname is constructed as primaryAccession.start.stop to match the SILVA FASTA sequence ID.

Usage:
  python build_taxtastic_seqinfo.py taxmap.txt seq_info.csv tax_ids.txt [--seqids file.txt]
"""

import sys
import csv
import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Build seq_info.csv and tax_ids.txt from SILVA taxmap (NCBI taxonomy IDs)'
    )
    parser.add_argument('taxmap', help='SILVA taxmap file (taxmap_slv_ssu_ref_nr_138.2.txt)')
    parser.add_argument('seq_info', help='Output seq_info.csv for taxit create')
    parser.add_argument('tax_ids', help='Output tax_ids.txt (unique NCBI tax_ids, one per line)')
    parser.add_argument('--seqids', metavar='FILE',
                        help='Text file of sequence IDs to include (one per line); '
                             'if omitted, all sequences in taxmap are included')
    args = parser.parse_args()

    allowed = None
    if args.seqids:
        with open(args.seqids) as f:
            allowed = set(line.strip() for line in f if line.strip())
        print(f"Filtering to {len(allowed)} sequence IDs from {args.seqids}")

    seq_rows = []
    tax_id_set = set()
    skipped = 0

    print(f"Reading {args.taxmap}...")
    with open(args.taxmap) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for i, row in enumerate(reader):
            accession = (row.get('primaryAccession') or '').strip()
            start = (row.get('start') or '').strip()
            stop = (row.get('stop') or '').strip()
            organism = (row.get('organism_name') or '').strip()
            tax_id = (row.get('ncbi_taxid') or row.get('taxid') or '').strip()

            # Construct seqname from parts (matches SILVA FASTA ID format)
            if accession and start and stop:
                seqname = f"{accession}.{start}.{stop}"
            else:
                skipped += 1
                continue

            if not tax_id:
                skipped += 1
                continue

            if allowed is not None and seqname not in allowed:
                continue

            seq_rows.append((seqname, accession or seqname, tax_id, organism, 'FALSE'))
            tax_id_set.add(tax_id)

            if (i + 1) % 100000 == 0:
                print(f"  Processed {i + 1} rows...")

    print(f"  {len(seq_rows)} sequences mapped, {skipped} skipped (missing seqname or tax_id)")

    if allowed is not None:
        missing = len(allowed) - len(seq_rows)
        if missing > 0:
            print(f"  WARNING: {missing} seqids had no entry in taxmap")

    print(f"Writing {args.seq_info}...")
    with open(args.seq_info, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['seqname', 'accession', 'tax_id', 'species_name', 'is_type'])
        writer.writerows(seq_rows)

    print(f"Writing {args.tax_ids} ({len(tax_id_set)} unique NCBI tax_ids)...")
    with open(args.tax_ids, 'w') as f:
        for tid in sorted(tax_id_set, key=lambda x: int(x) if x.isdigit() else 0):
            f.write(tid + '\n')

    print(f"Done.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
