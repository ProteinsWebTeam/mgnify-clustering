#!/bin/bash

if [ $# -le 0 ]
        then echo "Path not given for references"
else
    FOLDERDB=$1
    mkdir -p $FOLDERDB
    cd $FOLDERDB
    
    echo "verify Uniprot release version"

    if [[ $(head -1 /ebi/ftp/pub/databases/uniprot/relnotes.txt) == $(head -1 relnotes.txt) ]]; 
    then
        echo "Uniprot hasn't been updated"
        exit
    fi

    echo "Updating Uniprotkb"
    #put together the content of the SwissProt and TrEMBL files to obtain uniprotkb
    UNIPROT_VERSION=`cut -d ' ' -f3 <(head -1 /ebi/ftp/pub/databases/uniprot/relnotes.txt)`

    zcat /ebi/ftp/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz /ebi/ftp/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz > "uniprotkb_${UNIPROT_VERSION}.fasta" 
    cp /ebi/ftp/pub/databases/uniprot/relnotes.txt .

    #save uniprot accessions in separate file
    #cut -d'|' -f2 <(grep '>' "uniprotkb_${UNIPROT_VERSION}.fasta") > "uniprot_acc_${UNIPROT_VERSION}.txt"
    
fi

