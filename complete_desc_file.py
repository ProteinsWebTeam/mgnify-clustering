#!/usr/bin/env python3

###################
#
# @author T. Paysan-Lafosse
#
# @brief This script verifies if the pfbuild step has succedeed and runs extra scripts to complete the DESC file
#       It should be run from the Pfam directory containing the files
#
###################

import time
import os
import re


def check_pfambuild():
    if os.path.isfile("PFAMOUT") == True and os.path.getsize("PFAMOUT") != 0:
        return True

    elif os.path.isfile("pfbuild.log") and os.path.getsize("pfbuild.log") != 0:
        with open("pfbuild.log") as f:
            datafile = f.readlines()
            # verify that the process has finished running
            if "Resource usage summary:\n" in datafile:
                r = re.compile(
                    r"cannot create temp file for here-document: No space left on device"
                )
                newlines = list(filter(r.search, datafile))
                if len(newlines) > 0:  # if error found in the process
                    print(
                        "An error occured when running pfbuild, run it again, please check the log file (pfbuild.log)"
                    )
                    print(newlines)
                    return

    return False


def complete_desc():

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
    os.system(f"/nfs/production/xfam/pfam/software/bin/pqc-overlap-rdb.pl .")
    print("--- Completed in %.2f minutes ---" % ((time.time() - start_time) / 60))

    os.system(f"chmod -R g+w *")


if __name__ == "__main__":
    print("Verify if Pfam build has completed successfully")
    done = check_pfambuild()
    if done == True:
        complete_desc()
    elif done == False:
        print("Pfbuild hasn't completed yet")
