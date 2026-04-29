# Building a SILVA 138.2 Phylogenetic Tree

A standalone Python implementation of Ben Kaehler's workflow for building phylogenetic trees from SILVA 138.2 NR99 reference sequences.

## Quick Start

### 1. Install dependencies

```bash
# Install Python packages
pip install -r requirements.txt

# Install FastTree (required for tree building)
# macOS (uses faster VeryFastTree):
brew install veryfasttree

# Ubuntu/Debian:
sudo apt-get install fasttree

# CentOS/RHEL:
sudo yum install fasttree
```

### 2. Download SILVA 138.2 alignment

The SILVA 138.2 SSURef NR99 **ALIGNED** alignment (~1.4 GB gzipped, ~25 GB uncompressed):

```bash
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz
```

This will create a file named `SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz`

⚠️ **IMPORTANT**: Make sure to download the `full_align` version (aligned sequences), NOT the unaligned sequences.

**Alternative**: If wget is not available, download the file manually from:
https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz

### 3. Run the script

```bash
# Basic usage (4 threads by default)
python build_silva_tree.py SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz

# With more threads for faster tree building
python build_silva_tree.py SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz -t 16

# Specify output directory
python build_silva_tree.py SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz -o output/

# Keep intermediate masked alignment file
python build_silva_tree.py SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz --keep-masked
```

## What the script does

1. **Load alignment** - Reads SILVA 138.2 SSURef NR99 sequences (510,495 sequences)
2. **Mask gappy columns** - Removes alignment columns with >99.56% gaps (~2,350 columns kept)
3. **Export masked alignment** - Saves filtered alignment to `silva-138.2-masked-aln.fasta`
4. **Build tree** - Uses FastTree to construct the phylogenetic tree
5. **Root tree** - Creates a midpoint-rooted version for proper phylogenetic interpretation

## Output files

- `silva-138.2-tree.nwk` - Unrooted phylogenetic tree (Newick format)
- `silva-138.2-rooted-tree.nwk` - Midpoint-rooted phylogenetic tree (Newick format)
- `silva-138.2-masked-aln.fasta` - Masked alignment (removed by default, use `--keep-masked` to preserve)

## Visualization

View your tree online or with desktop software:

- **iTOL** (online): https://itol.embl.de/ - Upload `silva-138.2-rooted-tree.nwk`
- **FigTree** (desktop): http://tree.bio.ed.ac.uk/software/figtree/
- **Dendroscope** (desktop): http://dendroscope.org/

For large trees (430K+ tips), consider:
- Subsampling tips before visualization
- Using iTOL's interactive features
- Network-based viewers like Archaeopteryx

## Performance notes

Tree building with FastTree on the full alignment takes:
- ~9 hours with 1 thread
- ~1-2 hours with 8-16 threads
- Depends on CPU and available memory

Memory usage:
- Alignment loading: ~5 GB
- Gap counting: ~2 GB
- Tree building: ~8-10 GB

## Troubleshooting

**"FastTree command not found"**
- Install FastTree using the commands above

**"Out of memory"**
- The script reads the full 1.4 GB alignment into memory
- Ensure you have at least 16 GB RAM available
- Process on a machine with sufficient memory

**"Too many sequences"**
- To create a smaller test tree, you can:
  1. Use `head` to extract first N sequences: `zcat SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz | head -20000 > test.fasta`
  2. Run the script on the test file
  3. Verify the workflow before processing the full dataset

## References

- Original notebook: [silva-trees by Ben Kaehler](https://github.com/BenKaehler/silva-trees)
- SILVA database: https://www.arb-silva.de/
- FastTree documentation: http://microbesonline.org/fasttree/
- scikit-bio: http://scikit-bio.org/

## Citation

If you use this workflow, please cite:
- SILVA database: Quast et al. (2012) - https://doi.org/10.1093/nar/gks1195
- FastTree: Price et al. (2009) - https://doi.org/10.1371/journal.pone.0009490
- scikit-bio: McDonald et al. (2021) - https://joss.theoj.org/papers/10.21105/joss.02738
