#!/usr/bin/env python3

###################
#
# @author T. Paysan-Lafosse
#
# @brief This script generates MSA for clusters not matching Pfam db
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
        for c in range(1, nb_zeros):
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

        # Check whether directory already exists. If not then get to work.
        if (
            not os.path.isdir(familydir)
            and not os.path.isdir(donefamily)
            and not os.path.isdir(ignorefamily)
        ):
            # build good quality pfam family
            build_families.build_family(cluster_file, cluster_align, cluster_align_dir)
        else:
            print(f"Family {cluster_align} ({cluster_rep}) ignored or already processed")
            print(os.path.isdir(familydir))

        text = f"{cluster_align}\t{line}\n"
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

    args = parser.parse_args()
    # clusters_to_align_file=mgy_seqs.cluster_seq.fa_percent_mgnify_2+_no_pfam_sorted

    al = alignments()
    al.datadir = os.path.dirname(os.path.realpath(args.inputfile))
    # al.aligned_dir = os.path.join(al.datadir,"Pfam-M")
    al.aligned_dir = os.path.join(al.datadir, "Pfam-M_test")
    al.cluster_dir = os.path.join(al.datadir, "clusters")

    pfam_m_names = os.path.join(al.datadir, "corresponding_clusters.txt")

    # Load the pfam configuration
    al.load_pfam_config()

    print("Deleting old alignment files")
    try:
        os.remove(pfam_m_names)
    except FileNotFoundError:
        pass

    if args.deletedata != "No":
        print(f"Deleting old data {al.aligned_dir}")
        shutil.rmtree(al.aligned_dir)

    os.makedirs(al.aligned_dir, exist_ok=True)

    print("Starting clusters alignments")
    count = 1

    with open(args.inputfile, "r") as f, open(pfam_m_names, "w") as pf:
        pf.write("pf_id\tcluster_rep\tnb_seq\tpercent_mgnify\tperecent_swissprot\n")
        for line in f:
            if count < 2:
                pf.write(al.get_alignment(count, line))
                count += 1
            else:
                break

    print(f"Pfam building done, {count-1} clusters saved")

