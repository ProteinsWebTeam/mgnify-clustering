#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script filters truncated sequences from a fasta file

@arguments [-f INPUT_FILE]: fasta file with sequences to remove
           
"""


import argparse
import os
import re

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input_file", help="input file", required=True)
    args = parser.parse_args()

    if os.path.isfile(args.input_file):
        with open(args.input_file, "r") as input, open(f"{args.input_file}_clear", "w") as output:

            line = input.readline()
            count = 0
            while line:
                if re.search(r"^>MGY", line):
                    count = 0
                    if "PL=00" in line or "PL=01" in line or "PL=10" in line:
                        output.write(line)
                        count = 1
                elif count == 1:
                    output.write(line)
                line = input.readline()
