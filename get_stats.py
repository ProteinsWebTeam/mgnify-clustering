#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script split the clusters into multiple files 
        and counts the number of sequences per cluster and the percentage of MGnify and SwissProt sequences

@arguments [-i input_file]: file countaining clustered sequences (fasta)
            [-f protein_file]: file countaining the list of proteins not found in the Pfam database (can be empty or non existing)
            [-u username]: username for database connection
            [-p password]: password for database connection
            [-s schema]: schema for database connection (VIPREAD)
           
"""


import argparse
import os
import sys
import re
import cx_Oracle
import traceback


def getConnection(user, password, schema):
    """
    Set database connection
    """
    print("Connection to the database")
    connectString = "".join([user, "/", password, "@", schema])
    connection = cursor = None
    try:
        connection = cx_Oracle.connect(connectString)
        cursor = connection.cursor()
    except:
        stackTrace = traceback.format_exc()
        print(stackTrace)
        if "invalid username" in stackTrace:
            print("Could not connect to {0} as user {1}".format(schema, user))
            print(
                "NB if your oracle username contains the '$' character either escape it or surround it with quotes"
            )
            print('eg "ops$craigm" or ops\$craigm')
            print("Otherwise the shell will remove the '$' and all subsequent characters!")
            sys.exit(1)
        else:
            print(stackTrace)
            sys.exit(1)

    return connection, cursor


def get_proteins_not_in_pfam(user, password, schema, protein_file):
    connection, cursor = getConnection(user, password, schema)

    print("Searching for proteins not matching Pfam in the database")
    sql = """SELECT PROTEIN_AC FROM INTERPRO.PROTEIN
            MINUS
            SELECT PROTEIN_AC FROM INTERPRO.MATCH PARTITION (MATCH_DBCODE_H)
        """

    cursor.execute(sql)
    protein_list = [row[0] for row in cursor]
    # print(protein_list)
    with open(protein_file, "w") as f:
        for protein in protein_list:
            f.write(f"{str(protein)}\n")

    cursor.close()
    connection.close()

    return protein_list


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input_file", help="file countaining cluster_sequences (fasta)", required=True
    )
    parser.add_argument(
        "-f",
        "--protein_file",
        help="file countaining protein accessions not found in Pfam",
        required=True,
    )
    parser.add_argument("-u", "--user", help="username for database connection", required=True)
    parser.add_argument("-p", "--password", help="password for database connection", required=True)
    parser.add_argument("-s", "--schema", help="database schema to connect to", required=True)
    args = parser.parse_args()

    # get list of UniProt accessions not found in Pfam
    if not os.path.isfile(args.protein_file) or os.path.getsize(args.protein_file) == 0:
        protein_list = get_proteins_not_in_pfam(
            args.user, args.password, args.schema, args.protein_file
        )
    else:
        print("Loading proteins accessions into memory")
        with open(args.protein_file, "r") as f:
            protein_list = [line.strip("\n") for line in f]

    print("Getting clusters' statistics")
    if os.path.isfile(args.input_file):
        dirname = os.path.dirname(args.input_file)
        with open(args.input_file, "r") as input, open(
            f"{args.input_file}_percent_mgnify_all", "w"
        ) as output, open(f"{args.input_file}_percent_mgnify_2+_no_pfam", "w") as output2:
            line = input.readline()
            count = 0
            countmgy = 0
            countswiss = 0
            rep = ""
            content = ""
            clusterdir = os.path.join(dirname, "clusters")
            os.makedirs(clusterdir, exist_ok=True)
            nb_clusters = 0

            while line and nb_clusters < 10000:
                line = line.strip("\n")
                if re.search(r"^>[A-Z0-9]+$", line):
                    # print(line)
                    # saving data count in file when finding a new cluster
                    if rep != "":
                        # print(rep)
                        # generate % MGnify sequences found in cluster
                        percentmgy = round(countmgy * 100 / count, 2)
                        # generate % SwissProt sequences found in cluster
                        percentswiss = round(countswiss * 100 / count, 2)

                        output.write(f"{rep}\t{count}\t{percentmgy}\t{percentswiss}\n")

                        # create sub directories to save clusters files if more that 1 seq per cluster
                        # and at least 1 UniProt seq and 1 mgnify seq
                        if count > 1 and percentmgy < 100.0 and percentmgy > 0.0:
                            inpfam = False
                            if rep[0:3] == "MGY":
                                subdir = os.path.join(clusterdir, rep[0:8])
                            else:
                                subdir = os.path.join(clusterdir, rep[0:3])
                                if rep in protein_list:
                                    inpfam = True
                            # save if representative seq not found in pfam
                            if inpfam == False:
                                print(rep)
                                nb_clusters += 1
                                output2.write(f"{rep}\t{count}\t{percentmgy}\t{percentswiss}\n")
                                os.makedirs(subdir, exist_ok=True)
                                # save clusters in different files
                                with open(os.path.join(subdir, f"{rep}.fa"), "w") as clusterf:
                                    clusterf.write(content)

                    # reinitialise counters for next cluster
                    count = 0
                    countmgy = 0
                    countswiss = 0
                    content = ""
                    rep = line[1:]

                else:
                    if re.match(r">", line):  # line indicating new cluster
                        if re.match(r">sp", line):  # seq from SwissProt
                            countswiss += 1
                        elif re.match(r">MGY", line):  # seq from MGnify
                            countmgy += 1
                        # total number of seq
                        count += 1

                    # save line content to write in cluster files
                    content += f"{line}\n"

                line = input.readline()

            # compute and save last cluster
            print(f"nb_clusters: {nb_clusters}")
            if count != 0:
                percentmgy = round(countmgy * 100 / count, 2)
                percentswiss = round(countswiss * 100 / count, 2)

                output.write(f"{rep}\t{count}\t{percentmgy}\t{percentswiss}\n")

                if count > 1 and percentmgy < 100 and percentmgy > 0:
                    inpfam = False
                    if rep[0:3] == "MGY":
                        subdir = os.path.join(clusterdir, rep[0:8])
                    else:
                        subdir = os.path.join(clusterdir, rep[0:3])
                        if rep in protein_list:
                            inpfam = True
                    # save if representative seq not found in pfam
                    if inpfam == False:
                        output2.write(f"{rep}\t{count}\t{percentmgy}\t{percentswiss}\n")
                        os.makedirs(subdir, exist_ok=True)
                        with open(os.path.join(subdir, f"{rep}.fa"), "w") as clusterf:
                            clusterf.write(content)
            else:
                sys.exit(0)
