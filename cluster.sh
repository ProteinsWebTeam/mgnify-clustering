#!/usr/local/bin/bash

# @author T. Paysan-Lafosse
# @brief this script concatains fasta files and generates MGnify/UniProt-KB clusters
#usage: bsub -q production-rh74 -M 600000 -R "rusage[mem=600000]" -oo clustering_full_bidir.log -J cluster_full_bidir -Pbigmem -n 16 ./cluster.sh clustering.cfg

if [[ $# -lt 1 ]]; then
    echo "Illegal number of parameters. Usage: cluster.sh config_file"
    exit 2
fi

CONFIG_FILE=$1
if [[ -s $CONFIG_FILE ]];then
    . $CONFIG_FILE
else
    echo "Config file ${CONFIG_FILE} not found"
    exit
fi

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

set -e

DATADIR="${SCRIPTDIR}/data"
mkdir -p $DATADIR

#get latest version of UniProtKB
if [[ $UPDATE_UNIPROT == "yes" ]];then
    echo "Updating UniProt"
    bsub -J "get_uniprotkb" -oo "${DATADIR}/uniprot.log" "${SCRIPTDIR}/update_uniprotkb.sh" $DATADIR
fi

cd $DATADIR

#generate fasta file with MGnify and UniProt data
#MGNIFYDIR="/nfs/production/interpro/metagenomics/peptide_db/releases"
MGNIFY_REL="${MGNIFYDIR}/${MGNIFY_VERSION}"
MGNIFY_FASTA="${DATADIR}/${MGNIFY_VERSION}.fa"
MGNIFY_FASTA_FULL="${DATADIR}/${MGNIFY_VERSION}_clear.fa"

if [[ -d $MGNIFY_REL ]]; then
    if [[ -s $MGNIFY_FASTA_FULL ]]; then
        echo "${MGNIFY_FASTA_FULL} already exists"
    else
        echo "Concataining MGNIFY data"
        ls -d ${MGNIFY_REL}/mgy_proteins_* | xargs zcat > $MGNIFY_FASTA
        echo "Filtering full length sequences"
        python3 "${SCRIPTDIR}/filter_partial.py" -f $MGNIFY_FASTA
        mv "${DATADIR}/${MGNIFY_VERSION}.fa_clear" $MGNIFY_FASTA_FULL
        rm $MGNIFY_FASTA
    fi
else
    echo "${MGNIFY_REL} not found"
    exit
fi

if [[ $UPDATE_UNIPROT == "yes" ]];then
    bwait -w "ended(get_uniprotkb)"
fi

#concat MGnify and UniProtKB files
UNIPROT_VERSION=`cut -d ' ' -f3 <(head -1 "${DATADIR}/relnotes.txt")`
UNIPROT_FASTA="${DATADIR}/uniprotkb_${UNIPROT_VERSION}.fasta"
MGYUNI_FASTA="${DATADIR}/mgy_${MGNIFY_VERSION}_uniprot_${UNIPROT_VERSION}_full.fa"

if [[ -s $MGYUNI_FASTA ]]; then
    echo "${MGYUNI_FASTA} already exists"
else
    echo "Cat Mgnify and UniProt files"
    `cat $UNIPROT_FASTA $MGNIFY_FASTA_FULL > $MGYUNI_FASTA `  
    rm $UNIPROT_FASTA #$MGNIFY_FASTA_FULL 
fi

#Clustering
mmseqs="/nfs/production/interpro/metagenomics/peptide_db/mmseqs"
new_mmseqs="${SCRIPTDIR}/MMseqs2/build/bin/mmseqs"

cd $PFAM_DIR
SUBDIR="${MGNIFY_VERSION}_${UNIPROT_VERSION}_FULL_bidir"
mkdir -p $SUBDIR

DB="${SUBDIR}/mgy_seqs"
cluster_file="${DB}.cluster_seq.fa"

if [[ ! -s $cluster_file ]];then
    echo "Starting clustering"
    $mmseqs createdb $MGYUNI_FASTA $DB.mmseqs
    $mmseqs linclust $DB.mmseqs $DB.cluster "${SUBDIR}/tmp/" --min-seq-id 0.5 -c 0.75 --cov-mode 0 --threads 16
    $mmseqs result2repseq $DB.mmseqs $DB.cluster $DB.cluster_rep
    $mmseqs result2flat $DB.mmseqs $DB.mmseqs $DB.cluster_rep $DB.cluster_rep.fa
    $mmseqs createtsv $DB.mmseqs $DB.mmseqs $DB.cluster $DB.cluster.tsv
    $mmseqs createseqfiledb $DB.mmseqs $DB.cluster $DB.cluster_seq
    $mmseqs result2flat $DB.mmseqs $DB.mmseqs $DB.cluster_seq $DB.cluster_seq.fa
    echo "Clustering completed successfully"
fi

list_accessions="${cluster_file}_percent_mgnify_2+_no_pfam"

if [[ -s $cluster_file ]];then
    #counting cluster size and Mgnify percentages per cluster
    source ~oracle/ora112setup.sh
    prot_not_in_pfam="${DATADIR}/proteins_not_in_pfam_${UNIPROT_VERSION}.txt"
    if [[ -d "${SUBDIR}/clusters" ]]
    then
        echo "Deleting old clusters files"
        rm -r "${SUBDIR}/clusters"
    fi
    if [[ -s $list_accessions ]]
    then
        rm $list_accessions
    fi
    python "${SCRIPTDIR}/get_stats.py" -i $cluster_file -f $prot_not_in_pfam -u $USERNAME -p $PASSWORD -s $SCHEMA
else
    echo "Clustering failed"
    exit
fi

list_accessions_sorted="${list_accessions}_sorted"
if [[ -s $list_accessions ]];then
    sort -r -k2 -n $list_accessions > $list_accessions_sorted
else
    echo "Getting statistics failed or no cluster found with more than 2 sequences and containing a UniProt sequence not in Pfam"
fi


