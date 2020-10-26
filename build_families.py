#!/usr/bin/env python3

###################
#
# Script to build new families from Mnigfy-UniProt clutering
#
# Consider adding in code to help with curation such as
# Checking if any structures exist
# Running duffem.pl
# Add SE line to DESC
#
###################

from complete_desc_file import check_pfambuild, complete_desc
import os
from shutil import copyfile, rmtree
import argparse
import time
import re
import subprocess
import sys


def check_lift_over(family, outfile, count):

    family_seed = f"{family}_SEED"

    if os.path.isfile(outfile) and os.path.getsize(outfile) != 0:
        copyfile(outfile, "SEED4")
        print(f"Liftover complete for {family}")
        prepare_seed()

        if os.path.isfile("SEED") and os.path.getsize("SEED") != 0:
            print(f"Building Pfam from seed alignment for {family}")
            # os.system("pfbuild -withpfmake SEED")  # ~3 minutes run
            return True
        else:
            print(f"An error occured while building SEED alignment for {family}")
            if count < 1:
                return
            else:
                return "failed"

    elif os.path.isfile("liftover.log") and os.path.getsize("liftover.log") != 0:
        with open("liftover.log") as f:
            datafile = f.readlines()
            # verify that the process has finished running
            if "Resource usage summary:\n" in datafile:
                r = re.compile("Exited with exit code")
                newlines = list(filter(r.match, datafile))
                if len(newlines) > 0:  # if error found in the process
                    print(newlines)
                    if count < 1:  # number of failures below 1, try to run a second time
                        to_delete = ["*.fa", "*.log", "*.aln", "*.hmm", "hmmsearch.tbl"]
                        for filed in to_delete:
                            try:
                                os.remove(filed)
                            except FileNotFoundError:
                                pass
                        # Memory limit issue, stop processing this family
                        if "Exited with exit code 25.\n" in newlines:
                            print("Memory limit error")
                            return
                        else:
                            print("Error in liftover, trying again")
                            os.system(
                                f"perl /nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/liftover_alignment.pl -align {family_seed}"
                            )
                            return
                    else:
                        print(newlines)
                        return "failed"

    return False


def prepare_seed():
    trim_alignment_script = (
        "/nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/trim_alignment.pl"
    )
    # Make non redundant at 80% identity
    os.system("belvu -n 80 -o mul SEED4 | grep -v '//' > SEED3")

    # Remove gappy columns from N and C terminus
    print("Starting trim alignment")
    start_time = time.time()
    os.system(f"perl {trim_alignment_script} -in SEED3 -out SEED2")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    # Remove partial sequences, i.e. those with a gap character at start of end of alignment
    os.system("belvu -P -o mul SEED2 | grep -v '//' > SEED")


def build_alignment(cluster_file, family):

    liftover_alignment = (
        "/nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/liftover_alignment.pl"
    )
    alignment_script = (
        "/nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/create_alignment.pl"
    )

    print(f"Building family {family}")

    # align sequences in the cluster
    start_time = time.time()
    print(f"Starting alignment for {family}")
    if os.path.isfile(cluster_file):
        os.system(f"perl {alignment_script} -fasta {cluster_file} -m > {family}")
        print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

        # transform from aligned fasta to Pfam alignment format
        family_seed = f"{family}_SEED"
        while os.path.isfile(family_seed) == False or os.path.getsize(family_seed) == 0:
            os.system(f"belvu -o mul {family} | grep -v '//' >  {family_seed}")

        # align against pfamseq database
        print(f"Starting Liftover")
        os.system(f"perl {liftover_alignment} -align {family_seed}")  # ~7 minutes run

    else:
        print(f"File not found error: {cluster_file}")
        sys.exit()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--inputfile",
        help="path to fasta file containing the cluster sequences",
        required=True,
    )
    parser.add_argument("-f", "--family", help="Family accession to process", required=True)
    parser.add_argument("-d", "--datadir", help="Directory where to save data", required=True)

    args = parser.parse_args()

    scriptdir = os.path.dirname(os.path.realpath(__file__))
    familydir = os.path.join(args.datadir, args.family)
    outfile = os.path.join(familydir, f"{args.family}_SEED.phmmer")

    try:
        rmtree(familydir)
    except FileNotFoundError:
        pass

    os.makedirs(familydir, exist_ok=True)
    os.chdir(familydir)

    build_alignment(args.inputfile, args.family)

    # wait for lift_over to complete
    count_failed = 0
    success_liftover = False
    while True:
        done = check_lift_over(args.family, outfile, count_failed)
        if done == True:  # liftover completed successfully
            success_liftover = True
            break
        elif done == "failed":  # tried but failed running liftover, trying again
            count_failed += 1
        elif done == None:  # tried but failed running liftover with Memory limit
            count_failed += 1
        else:
            continue

        if count_failed > 1:
            break
        time.sleep(0.2)

    # if lift_over successfully completed, carry on with pfbuild
    # success_liftover = True
    if success_liftover:
        pf = subprocess.Popen(
            ["pfbuild", "-withpfmake", "SEED"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        err = pf.communicate()[1].decode("utf-8")
        if pf.returncode != 0:
            print(f"Error while running pfbuild: {err}")
            sys.exit()

        count_failed = 0
        # wait for pfbuild to complete
        while True:
            done = check_pfambuild()
            if done == True:  # pfbuild completed, complete DESC file
                complete_desc()
                break
            elif done == None:  # pfbuild failed to complete
                if count_failed < 1:
                    os.remove("pfbuild.log")
                    os.system("pfbuild -withpfmake SEED")
                count_failed += 1
            else:
                continue

            if count_failed > 1:
                break

            time.sleep(0.2)

    os.chdir(scriptdir)
    os.system(f"chmod -R g+w {familydir}")
