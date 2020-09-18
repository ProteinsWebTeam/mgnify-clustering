# MGnify clustering
Pipeline finding Mgnify/UniProt clusters to increase the metagenomic coverage in Pfam

## Requirements
Python 3.6 and above with Cx_Oracle package
MMseqs2 (https://github.com/soedinglab/MMseqs2)

Config file with the variables as specified in `clustering_model.cfg`

## Clustering
The pipeline is divided in 2 steps, both executed from `cluster.sh`:
- Clustering 
- Getting clusters with more than 2 sequences including MGnify and UniProt sequences where the representative sequence isn't found in Pfam

Usage: `bsub -q production-rh74 -M 600000 -R "rusage[mem=600000]" -oo clustering_full_bidir.log -J cluster_full_uni -Pbigmem -n 16 ./cluster.sh clustering.cfg`

## Buiding family
This is the next step after the clustering, allowing the automatic build of good quality Pfam families.
It runs a series of perl scripts in order to get a good quality family. It needs to be executed on an interactive shell, incompatible with `bsub` command.

- Generation of the multiple sequence alignment for the clusters selected in the previous step
- Alignment against pfamseq database
- Keeping only non redundant at 80% identity sequences
- Removing gappy columns from N and C terminus
- Removing partial sequences, i.e. those with a gap character at start of end of alignment
- Running pfbuild to generate missing files
- Searching for PDB references
- Searching for SwissProt information
- Searching for matching species
- Adding DUF identifier
- Adding family identifier
- Running final checks

Usage: `python generate_alignments.py -i file_containing_stats_sorted_by_cluster_size_reverse [-d yes (delete previous files)] [-b family to start building from] [-n number of families to build]`