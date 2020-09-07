#!/usr/bin/env python3

###################
#
# @author T. Paysan-Lafosse
#
# @brief This script generates MSA for clusters not matching Pfam db
#
# usage: python generate_alignments.py -i file_containing_stats_sorted_by_cluster_size_reverse 
#                                     [-d "yes" (delete previous files)] 
#                                     [-b family to start building from] 
#                                     [-n number of families to build]
#
###################

import os
import sys
import shutil
import argparse
import build_families


class alignments:
    def __init__(self):
        self.cluster_dir = ""
        self.aligned_dir = ""
        self.datadir = ""

    def chunks(self, lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    def load_pfam_config(self):
        if os.path.isfile("/nfs/production/xfam/pfam/pfamrc"):
            os.system("source /nfs/production/xfam/pfam/pfamrc")
        else:
            print("Can't find pfamrc file")
            sys.exit()

    def get_cluster_file(self, cluster_rep):
        if cluster_rep[0:3] == "MGY":
            return os.path.join(self.cluster_dir, f"{cluster_rep[0:8]}/{cluster_rep}.fa")
        else:
            return os.path.join(self.cluster_dir, f"{cluster_rep[0:3]}/{cluster_rep}.fa")

    def get_cluster_align(self, count):
        nb_zeros = 6 - len(str(count))
        cluster_align = "Pfam-M_"
        for c in range(0, nb_zeros):
            cluster_align = f"{cluster_align}0"

        return f"{cluster_align}{count}"

    def get_alignment(self, count, line):
        line = line.strip()
        cluster_rep = line.split("\t")[0]

        cluster_file = self.get_cluster_file(cluster_rep)
        cluster_align = self.get_cluster_align(count)

        cluster_align_dir = os.path.join(self.aligned_dir, cluster_align)
        print(cluster_rep, cluster_align)

        familydir = os.path.join(self.aligned_dir, cluster_align)
        donefamily = os.path.join(self.aligned_dir, f"DONE/{cluster_align}")
        ignorefamily = os.path.join(self.aligned_dir, f"IGNORE/{cluster_align}")

        success = False
        # Check whether directory already exists. If not then get to work.
        if (
            not os.path.isdir(familydir)
            and not os.path.isdir(donefamily)
            and not os.path.isdir(ignorefamily)
        ):
            while success == False:
                # build good quality pfam family
                success = build_families.build_family(
                    cluster_file, cluster_align, cluster_align_dir
                )
            text = f"{cluster_align}\t{line}\n"
        else:
            print(f"Family {cluster_align} ({cluster_rep}) ignored or already processed")
            text = None
        print(text)
        return text


# verify input parameters are given
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputfile", help="file countaining cluster statistics", required=True
    )
    parser.add_argument(
        "-d",
        "--deletedata",
        help="Specify if you want to delete previously generated data (default=No)",
        default="No",
    )
    parser.add_argument(
        "-b",
        "--begin_count",
        help="Family number value to start the alignment from (default=1)",
        default=1,
    )
    parser.add_argument(
        "-n",
        "--number_to_process",
        help="Number of clusters to process (default=1000)",
        default=1000,
    )

    args = parser.parse_args()
    # clusters_to_align_file=mgy_seqs.cluster_seq.fa_percent_mgnify_2+_no_pfam_sorted

    al = alignments()
    al.datadir = os.path.dirname(os.path.realpath(args.inputfile))
    al.aligned_dir = os.path.join(al.datadir, "Pfam-M")
    # al.aligned_dir = os.path.join(al.datadir, "Pfam-M_test")
    al.cluster_dir = os.path.join(al.datadir, "clusters")

    pfam_m_names = os.path.join(al.datadir, "corresponding_clusters.txt")
    # pfam_m_names = os.path.join(al.datadir, "corresponding_clusters_test.txt")

    # Load the pfam configuration
    al.load_pfam_config()

    if args.deletedata != "No":
        print("Deleting old alignment files")
        try:
            os.remove(pfam_m_names)
        except FileNotFoundError:
            pass
        print(f"Deleting old data {al.aligned_dir}")
        shutil.rmtree(al.aligned_dir)

    print(al.aligned_dir)
    os.makedirs(al.aligned_dir, exist_ok=True)

    print("Starting clusters alignments")
    count = int(args.begin_count)
    if int(args.number_to_process) > 1:
        final_cluster = int(args.begin_count) + int(args.number_to_process)
    else:
        final_cluster = int(args.begin_count)
    print(f"Processing clusters {count} to {final_cluster}")

    with open(args.inputfile, "r") as f, open(pfam_m_names, "a") as pf:
        if args.deletedata != "No":
            pf.write("pf_id\tcluster_rep\tnb_seq\tpercent_mgnify\tperecent_swissprot\n")

        tmp_count = 0  # used if do not want to start building from first cluster
        for line in f:
            tmp_count += 1
            if tmp_count >= count:  # start building from given file line
                if count <= final_cluster:
                    text = al.get_alignment(count, line)
                    if text != None:
                        pf.write(text)
                    count += 1
                else:
                    break
            else:
                pass

    print(f"Pfam building done, {args.number_to_process} clusters saved")

