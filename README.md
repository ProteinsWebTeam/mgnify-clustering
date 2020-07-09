# MGnify clustering
Pipeline finding Mgnify/UniProt clusters to increase the metagenomic coverage in Pfam

## Requirements
Python 3.6 and above
Config file with the variables as specified in `clustering_model.cfg`

## Clustering
The pipeline is divided in 3 steps, all executed from `cluster.sh`:
- Clustering 
- Getting clusters with more than 2 sequences including MGnify and UniProt sequences where the representative sequence isn't found in Pfam
- Generating the multiple sequence alignment for the clusters selected in the previous step

Usage: `bsub -q production-rh74 -M 600000 -R "rusage[mem=600000]" -o clustering_full_bidir.log -J cluster_full_uni -Pbigmem -n 16 ./cluster.sh clustering.cfg`