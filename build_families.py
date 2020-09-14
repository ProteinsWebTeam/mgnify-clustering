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

import os
from shutil import copyfile, rmtree
import argparse
import time
import re


def build_family(cluster_file, family, familydir):
    root = "/homes/agb/Work/DB/PfamBMM/Pfam-B_33.1"
    liftover_alignment = (
        "/nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/liftover_alignment.pl"
    )
    alignment_script = (
        "/nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/create_alignment.pl"
    )
    trim_alignment_script = (
        "/nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/trim_alignment.pl"
    )

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    print(f"Building family {family}")

    try:
        rmtree(familydir)
    except FileNotFoundError:
        pass
    os.makedirs(familydir, exist_ok=True)
    os.chdir(familydir)

    # align sequences in the cluster
    start_time = time.time()
    print(f"Starting alignment for {family}")
    os.system(f"perl {alignment_script} -fasta {cluster_file} -m > {family}")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    # transform from aligned fasta to Pfam alignment format
    family_seed = f"{family}_SEED"
    while os.path.isfile(family_seed) == False or os.path.getsize(family_seed) == 0:
        os.system(f"belvu -o mul {family} | grep -v '//' >  {family_seed}")

    # align against pfamseq database
    print(f"Starting Liftover")
    outfile = f"{family_seed}.phmmer"
    success_liftover = False

    while success_liftover != True:
        os.system(f"perl {liftover_alignment} -align {family_seed}")  # ~7 minutes run
        while os.path.isfile(outfile) == False or os.path.getsize(outfile) == 0:
            time.sleep(30)
            if os.path.isfile("liftover.log") and os.path.getsize("liftover.log") != 0:
                with open("liftover.log") as f:
                    datafile = f.readlines()
                    # verify that the process has finished running
                    if "Resource usage summary:\n" in datafile:
                        r = re.compile("Exited with exit code")
                        newlines = list(filter(r.match, datafile))
                        if len(newlines) > 0:  # if error found in the process
                            print(newlines)
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
                                break
                    else:
                        success_liftover = True
        success_liftover = True
        print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    copyfile(outfile, "SEED4")
    print(f"Liftover complete")

    # Make non redundant at 80% identity
    os.system("belvu -n 80 -o mul SEED4 | grep -v '//' > SEED3")

    # Remove gappy columns from N and C terminus
    print("Starting trim alignment")
    start_time = time.time()
    os.system(f"perl {trim_alignment_script} -in SEED3 -out SEED2")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    # Remove partial sequences, i.e. those with a gap character at start of end of alignment
    os.system("belvu -P -o mul SEED2 | grep -v '//' > SEED")

    print("Building Pfam from seed alignment")
    success_pfbuild = False
    while success_pfbuild == False:
        start_time = time.time()
        os.system("pfbuild -withpfmake SEED")  # ~3 minutes run
        while os.path.isfile("PFAMOUT") == False or os.path.getsize("PFAMOUT") == 0:
            time.sleep(30)
            if os.path.isfile("pfbuild.log") and os.path.getsize("pfbuild.log") != 0:
                with open("pfbuild.log") as f:
                    datafile = f.readlines()
                    # verify that the process has finished running
                    if "Resource usage summary:\n" in datafile:
                        r = re.compile("Exited with exit code")
                        newlines = list(filter(r.match, datafile))
                        if len(newlines) > 0:  # if error found in the process
                            os.remove("pfbuild.log")
                            print("Error in pfbuild")
                            print(newlines)
                            break
                    else:
                        success_pfbuild = True
        success_pfbuild = True
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    print("Searching for PDB reference")
    start_time = time.time()
    os.system(f"perl /homes/agb/Scripts/add_pdb_ref.pl")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    print("Searching for SwissProt info")
    start_time = time.time()
    os.system(f"perl /nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/swissprot.pl -num 10")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    print("Searching for species")
    start_time = time.time()
    os.system(f"perl /homes/agb/Scripts/species_summary.pl .")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    print("Adding DUF identifier")
    start_time = time.time()
    os.system(f"perl /homes/agb/Scripts/duffem.pl -overwrite -duf .")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    print("Adding DUF extra step")
    start_time = time.time()
    os.system(f"perl /nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/nextDUF.pl")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    print("Running final checks")
    start_time = time.time()
    os.system(f"/nfs/production/xfam/pfam/software/bin/pqc-all .")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    os.chdir(scriptdir)
    os.system(f"chmod -R g+w {familydir}")

    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--family", help="Family accession to process", required=True)
    parser.add_argument(
        "-d", "--datadir", help="Directory containing the data to process", required=True
    )

    args = parser.parse_args()
    build_family(args.family, args.datadir)

