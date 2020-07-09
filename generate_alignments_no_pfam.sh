#!/usr/local/bin/bash

###################
#
# @author T. Paysan-Lafosse
#
# @brief This script generates MSA for clusters not matching Pfam db
#
###################

#pretty_sort () {
   # PrettySort means sorting by p = min(height, cons_len)*avg_cons

    #where height=nr of seqs, and cons_len and avg_cons come from piping belvu -c output to:

    # #!/usr/bin/gawk -f
    # if($4== "=") {
    #     n++;
    #     if ($7 > 0.5) {c++; }
    # }
    # (/Aver/) {a = $NF}
    # END {print "total_len=", n, "cons_len=", c, "avg_cons=", a}
#}

# verify input parameters are given
if [[ $# -lt 1 ]]; then
    echo "Illegal number of parameters. Usage: generate_alignments_no_pfam.sh cluster_list"
    exit 2
else
    clusters_to_align_file="$1"
fi

#Load the pfam configuration
if [ -f /nfs/production/xfam/pfam/pfamrc ]; then
  source /nfs/production/xfam/pfam/pfamrc
else
    echo "Can't find pfamrc file"
    exit
fi

alignment="/nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/create_alignment.pl"
#clusters_to_align_file=mgy_seqs.cluster_seq.fa_percent_mgnify_2+_no_pfam

dir=`dirname $clusters_to_align_file`
aligned_dir="${dir}/aligned"
mkdir -p $aligned_dir

while IFS= read -r line; do
    tmp_array=($(echo $line | tr "\t" "\n"))
    cluster_rep=${tmp_array[0]}

    if [[ ${cluster_rep:0:3} == "MGY" ]];then
        cluster_dir=${cluster_rep:0:8}
    else
        cluster_dir=${cluster_rep:0:3}
    fi
    cluster_file="${dir}/clusters/${cluster_dir}/${cluster_rep}.fa"

    mkdir -p "${aligned_dir}/${cluster_dir}"
    aligned_file="${aligned_dir}/${cluster_dir}/${cluster_rep}"
    
    perl $alignment -fasta $cluster_file -m > $aligned_file
# | xargs belvu -c | gawk -f
#     ($4== "=") {
#         n++;
#         if ($7 > 0.5) {c++; }
#     }
#     (/Aver/) {a = $NF}
#     END {print "total_len=", n, "cons_len=", c, "avg_cons=", a} `#| pretty_sort
done < $clusters_to_align_file