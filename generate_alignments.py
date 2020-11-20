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
import redis
import time
import subprocess
from multiprocessing import Process

import build_families
from complete_desc_file import complete_desc, check_pfambuild


class alignments:
    def __init__(self):
        self.cluster_dir = ""
        self.aligned_dir = ""
        self.datadir = ""
        self.scriptdir = os.path.dirname(os.path.realpath(__file__))
        self.queue_in = "lift_over_in"
        self.queue_done = "lift_over_done"
        self.pfam_in = "pfam_in"
        self.pfam_failed = dict()
        self.liftover_failed = dict()
        self.server = redis.StrictRedis(
            port=6379, host="localhost", db=0, charset="utf-8", decode_responses=True
        )

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

        familydir = os.path.join(self.aligned_dir, cluster_align)
        print(cluster_rep, cluster_align)

        donefamily = os.path.join(self.aligned_dir, f"DONE/{cluster_align}")
        donemfamily = os.path.join(self.aligned_dir, f"DONE_MERGED/{cluster_align}")
        ignorefamily = os.path.join(self.aligned_dir, f"IGNORE/{cluster_align}")
        functionfamily = os.path.join(self.aligned_dir, f"FUNCTION/{cluster_align}")
        failedfamily = os.path.join(self.aligned_dir, f"FAILED/{cluster_align}")
        duffamily = os.path.join(self.aligned_dir, f"DUF/{cluster_align}")

        # Check whether directory already exists. If not then get to work.
        if (
            not os.path.isdir(familydir)
            and not os.path.isdir(donefamily)
            and not os.path.isdir(ignorefamily)
            and not os.path.isdir(functionfamily)
            and not os.path.isdir(duffamily)
            and not os.path.isdir(donemfamily)
            and not os.path.isdir(failedfamily)
        ):
            text = f"{cluster_align}\t{line}\n"

            # build good quality pfam family
            os.makedirs(familydir, exist_ok=True)
            os.chdir(familydir)
            build_families.build_alignment(cluster_file, cluster_align)
            self.server.rpush(self.queue_in, cluster_align)
            os.chdir(self.scriptdir)
        else:
            print(
                f"Family {cluster_align} ({cluster_rep}) ignored, processing or already processed"
            )
            text = None
        return text

    def wait_liftover(self):
        print("waiting for liftover to complete")
        print(f"{self.server.llen(self.queue_in)} families to process")

        while True:
            family = self.server.lpop(self.queue_in)
            if not family:
                break
            # print(family)

            familydir = os.path.join(self.aligned_dir, family)
            os.chdir(familydir)
            faileddir = os.path.join(self.aligned_dir, "FAILED")

            outfile = os.path.join(familydir, f"{family}_SEED.phmmer")
            count = self.liftover_failed[family] if family in self.liftover_failed else 0
            done = build_families.check_lift_over(family, outfile, count)

            if done == False:  # liftover hasn't completed yet
                self.server.rpush(self.queue_in, family)
            elif done == True:  # liftover completed successfully
                pf = subprocess.Popen(
                    ["pfbuild", "-withpfmake", "SEED"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                err = pf.communicate()[1].decode("utf-8")
                if pf.returncode != 0:
                    print(f"Error while running pfbuild: {err}")
                    self.pfam_failed[family] = 1
                    os.system(f"mv {familydir} {faileddir}")
                else:
                    self.server.rpush(self.queue_done, family)
                    print(err)
            elif done == "failed":  # tried but failed running liftover twice
                self.liftover_failed[family] += 1
                os.system(f"mv {familydir} {faileddir}")
            else:  # tried but failed running liftover once, trying again
                self.liftover_failed[family] = 1
                self.server.rpush(self.queue_in, family)

            time.sleep(0.2)
        # print(f"All liftover completed" , {len(self.liftover_failed)} failed {self.liftover_failed.keys()}")

    def wait_pfbuild(self):
        print("waiting for pfbuild to complete")
        while True:
            if not self.server.llen(self.queue_in) and not self.server.llen(self.queue_done):
                break
            family = self.server.lpop(self.queue_done)
            if not family:
                time.sleep(0.2)
                continue
            # print(family)

            self.server.rpush(self.pfam_in, family)
            familydir = os.path.join(self.aligned_dir, family)
            faileddir = os.path.join(self.aligned_dir, "FAILED")
            os.chdir(familydir)
            done = check_pfambuild()

            if done == False:  # pfbuild hasn't completed yet
                self.server.rpush(self.queue_done, family)
            elif done == None:  # pfbuild failed to complete
                if family not in self.pfam_failed:  # attempt to run pfbuild a second time
                    os.remove("pfbuild.log")
                    os.system("pfbuild -withpfmake SEED")
                    self.server.rpush(self.queue_done, family)
                    self.pfam_failed[family] = 1
                else:
                    self.pfam_failed[family] += 1
                    os.chdir(self.scriptdir)
                    os.system(f"mv {familydir} {faileddir}")
            else:  # pfbuild completed successfully
                complete_desc()
                print(f"Family {family} successfully built")
                os.chdir(self.scriptdir)
                os.system(f"chmod -R g+w {familydir}")
            self.server.lpop(self.pfam_in)


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
        final_cluster = int(args.begin_count) + int(args.number_to_process) - 1
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

    liftover_failed = dict()
    liftover = Process(target=al.wait_liftover())
    liftover.start()

    pfam_failed = dict()
    pfambuild = Process(target=al.wait_pfbuild())
    pfambuild.start()

    liftover.join()
    pfambuild.join()

    os.chdir(al.scriptdir)

    logfile = "generate_alignments.log"

    with open(logfile, "a") as log:
        logtime = {time.strftime("%d %b %Y %H:%M", time.localtime())}
        log.write(f"\n{logtime}:\n")
        print(f"Pfam building done, {args.number_to_process} clusters processed")
        log.write(f"Pfam building done, {args.number_to_process} clusters processed\n")

        if len(al.liftover_failed) > 0:
            print(
                f"{len(al.liftover_failed)} failed liftover: {' '.join(al.liftover_failed.keys())}"
            )
            log.write(
                f"{len(al.liftover_failed)} failed liftover: {' '.join(al.liftover_failed.keys())}\n"
            )
        if len(al.pfam_failed) > 0:
            print(f"{len(al.pfam_failed)} failed pfbuild: {' '.join(al.pfam_failed.keys())}")
            log.write(f"{len(al.pfam_failed)} failed pfbuild: {' '.join(al.pfam_failed.keys())}\n")

