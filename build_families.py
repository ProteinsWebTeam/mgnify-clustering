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
from shutil import copyfile
import argparse
import time


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
    pdb_ref_script = "/homes/agb/Scripts/add_pdb_ref.pl"

    swissprot_script = "/nfs/production/xfam/pfam/software/Pfam/PfamScripts/make/swissprot.pl"

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    print(f"Building family {family}")

    os.makedirs(familydir, exist_ok=True)
    os.chdir(familydir)
    print(os.getcwd())

    # align sequences in the cluster
    print(f"Starting alignment for {family}")
    os.system(f"perl {alignment_script} -fasta {cluster_file} -m > {family}")

    # transform from aligned fasta to Pfam alignment format
    family_seed = f"{family}_SEED"
    os.system(f"belvu -o mul {family} | grep -v '//' >  {family_seed}")
    print(os.path.getsize(family_seed))

    # align against pfamseq database
    print(f"Starting Liftover")
    print(family_seed)
    os.system(f"perl {liftover_alignment} -align {family_seed}")
    outfile = f"{family_seed}.phmmer"

    while os.path.isfile(outfile) == False or os.path.getsize(outfile) == 0:
        time.sleep(30)

    copyfile(outfile, "SEED4")
    print(f"Liftover complete")

    # Make non redundant at 80% identity
    os.system("belvu -n 80 -o mul SEED4 | grep -v '//' > SEED3")

    # Remove gappy columns from N and C terminus
    print("Starting trim alignment")
    os.system(f"perl {trim_alignment_script} -in SEED3 -out SEED2")
    print("Trim alignment complete")

    # Remove partial sequences, i.e. those with a gap character at start of end of alignment
    os.system("belvu -P -o mul SEED2 | grep -v '//' > SEED")

    print("Building Pfam from seed alignment")
    os.system("pfbuild SEED")

    print("Searching for PDB reference")
    os.system(f"perl {pdb_ref_script}")

    print("Searching for SwissProt info")
    os.system(f"perl {swissprot_script}")

    os.chdir(scriptdir)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--family", help="Family accession to process", required=True)
    parser.add_argument(
        "-d", "--datadir", help="Directory containing the data to process", required=True
    )

    args = parser.parse_args()
    build_family(args.family, args.datadir)

