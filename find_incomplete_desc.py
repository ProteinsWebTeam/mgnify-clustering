#!/usr/bin/env python3

###################
#
# @author T. Paysan-Lafosse
#
# @brief This script checks and completes incomplete DESC files
#
# @args -d dir: directory containing the data to check
#
##################

import os
import subprocess
from complete_desc_file import check_pfambuild, complete_desc
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", help="directory containing the data to check", required=True)

    args = parser.parse_args()

    count = 0
    for family in os.listdir(args.dir):
        count += 1
        if count > 500 and count < 800:
            familydir = os.path.join(args.dir, family)
            desc_f = os.path.join(familydir, "DESC")
            if os.path.isfile(desc_f):
                proc = subprocess.Popen(["grep", "Who RU", desc_f], stdout=subprocess.PIPE)
                outs, errs = proc.communicate()
                outs = outs.decode("utf-8")
                if outs != "":
                    print(family)
                    os.chdir(familydir)
                    if check_pfambuild():
                        # print("pfbuild completed")
                        complete_desc()
                    else:
                        print("Pfbuild failed")
            else:
                print(f"No DESC file found for {family}")

