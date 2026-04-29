# SILVA SEPP-style phylogenetic tree for fragment insertion

## Author: Peter Falb

## First created April 29, 2026

The approach here is loosely based off of Ben Kaehler's approach to build a SILVA tree (<https://gist.github.com/BenKaehler/d9291d59bce5cd3d2a90c73b822b3a21>) The idea is to replicate the SEPP QIIME2 fragment insertion plugin, which used a tree constructed from the SILVA 128 NR99 database. Here we are updating the tree to the most recent version of SILVA, for use with an epa-ng and gappa based approach to fragment insertion.

## Basic approach:

1.  Imported SILVA SSURef NR99 full alignment (<https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz>)
2.  Convert from RNA to DNA using scikit bio
3.  Mask alignment columns that are very gappy (following Ben Kaehler's approach)
4.  Build a tree using fasttree

## Main Files

*build_silva_tree_clean.py* - Python file used to construct tree and mask alignments, etc

*submit_tree_job.sbatch* - SLURM script for submission to Agate supercomputer

*slurm_setup.sh* - setup script for virtual environment and module setup on Agate supercomputer

## Citations/Acknowledgements:

Basic approach replicates Ben Kaehler's approach on his above-linked Github, which itself attempts to replicate the approach taken for SEPP (<https://github.com/smirarab/sepp-refs/tree/master>)

*The SEPP Paper is here:*

Mirarab, S., Nguyen, N., & Warnow, T. (2011). SEPP: SATé-Enabled Phylogenetic Placement. In Biocomputing 2012 (pp. 247–258). WORLD SCIENTIFIC. <https://doi.org/10.1142/9789814366496_0024> 9789814366496_0024

*FastTree should also be cited:*

Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 5(3):e9490. <doi:10.1371/journal.pone.0009490>
