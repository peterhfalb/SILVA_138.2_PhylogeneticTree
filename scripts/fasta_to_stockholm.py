#!/usr/bin/env python3
"""
Convert a FASTA multiple sequence alignment to Stockholm format.

Each sequence is written as a single unbroken line. The output is compatible
with HMMER3 and Infernal tools, and suitable for use as the --aln-sto argument
to 'taxit create'.

Usage:
  python fasta_to_stockholm.py input.fasta output.sto
"""

import sys
import argparse


def convert(fasta_path, sto_path):
    seqs = {}
    order = []
    seqid = None

    print(f"Reading {fasta_path}...")
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                seqid = line[1:].split()[0]
                seqs[seqid] = []
                order.append(seqid)
            elif seqid:
                seqs[seqid].append(line)

    print(f"Writing {sto_path} ({len(order)} sequences)...")
    with open(sto_path, 'w') as f:
        f.write('# STOCKHOLM 1.0\n\n')
        for seqid in order:
            f.write(f'{seqid}  {"".join(seqs[seqid])}\n')
        f.write('//\n')

    print(f"Done.")
    return 0


def main():
    parser = argparse.ArgumentParser(description='Convert FASTA MSA to Stockholm format')
    parser.add_argument('fasta', help='Input FASTA alignment')
    parser.add_argument('stockholm', help='Output Stockholm (.sto) file')
    args = parser.parse_args()
    return convert(args.fasta, args.stockholm)


if __name__ == '__main__':
    sys.exit(main())
