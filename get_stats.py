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

from multiprocessing.dummy import Pool
import itertools


class process_cluster:
    def __init__(self):
        self.protein_dict = {}
        self.dirname = ""
        self.clusterdir = ""
        self.count = 0

    def getConnection(self, user, password, schema):
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

    def get_proteins_not_in_pfam(self, user, password, schema, protein_file):
        connection, cursor = self.getConnection(user, password, schema)

        print("Searching for proteins not matching Pfam in the database")
        sql = """SELECT PROTEIN_AC FROM INTERPRO.PROTEIN
                MINUS
                SELECT PROTEIN_AC FROM INTERPRO.MATCH PARTITION (MATCH_DBCODE_H)
            """

        cursor.execute(sql)
        self.protein_dict = [row[0] for row in cursor]
        # print(protein_dict)
        with open(protein_file, "w") as f:
            for protein in self.protein_dict:
                f.write(f"{str(protein)}\n")

        cursor.close()
        connection.close()

    def cluster_separator(self, line):
        if re.match(r"^>[A-Z0-9]+$", line):
            return line

    def process(self, chunk):
        inpfam = 0
        countswiss = 0
        countmgy = 0
        count_total = 0
        content = ""
        rep = ""
        counting = ""

        if self.count < 10000:  # end search if 10K clusters found
            for item in chunk:
                skip = False
                if re.match(r"^>", item):
                    if re.match(r">sp", item):  # seq from SwissProt
                        countswiss += 1
                        acc = item.split("|")[1]
                    elif re.match(r">MGY", item):  # seq from MGnify
                        countmgy += 1
                        acc = item.split(">")[1].split(" ")[0]
                        skip = True
                    else:  # seq from TrEMBL
                        acc = item.split("|")[1]

                    # print(acc)

                    if not skip and acc in self.protein_dict:
                        inpfam += 1
                        break

                    # total number of seq
                    count_total += 1

                    # save line content to write in cluster files
                    if count_total == 1:
                        content += f">{acc}\n"
                        rep = acc

                content += item

            # for clusters not found in Pfam and with at least 2 sequences
            if inpfam == 0 and count_total > 1:
                # print(f"Cluster {rep} not in Pfam")

                # generate % MGnify sequences found in cluster
                percentmgy = round(countmgy * 100 / count_total, 2)
                # generate % SwissProt sequences found in cluster
                percentswiss = round(countswiss * 100 / count_total, 2)

                # create sub directories to save clusters files if at least 1 UniProt seq and 1 mgnify seq
                if percentmgy < 100.0 and percentmgy > 0.0:
                    if rep[0:3] == "MGY":
                        subdir = os.path.join(self.clusterdir, rep[0:8])
                    else:
                        subdir = os.path.join(self.clusterdir, rep[0:3])
                    os.makedirs(subdir, exist_ok=True)

                    # save clusters in different files
                    with open(os.path.join(subdir, f"{rep}.fa"), "w") as clusterf:
                        clusterf.write(content)

                    counting = f"{rep}\t{count_total}\t{percentmgy}\t{percentswiss}\n"
                    self.count += 1
                    return counting


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

    pc = process_cluster()

    # get list of UniProt accessions not found in Pfam
    if not os.path.isfile(args.protein_file) or os.path.getsize(args.protein_file) == 0:
        pc.get_proteins_not_in_pfam(args.user, args.password, args.schema, args.protein_file)
    else:
        print("Loading proteins accessions into memory")
        with open(args.protein_file, "r") as f:
            pc.protein_dict = {line.strip("\n"): 1 for line in f}
    # print(f"protein_dict: {len(pc.protein_dict)}")

    pc.dirname = os.path.dirname(args.input_file)
    pc.clusterdir = os.path.join(pc.dirname, "clusters")
    os.makedirs(pc.clusterdir, exist_ok=True)

    # protein_dict = []
    print("Getting clusters' statistics")

    num_chunks = 5
    results = []

    if os.path.isfile(args.input_file):
        with open(args.input_file, "r") as inputf:
            chunks = itertools.groupby(inputf, pc.cluster_separator)
            while True:
                groups = []
                for key, chunk in itertools.islice(chunks, num_chunks):
                    if not key:
                        groups.append(list(chunk))
                if groups:
                    with Pool(8) as pool:
                        result = pool.map(pc.process, groups)
                        results.extend(result)
                else:
                    break

        with open(f"{args.input_file}_percent_mgnify_2+_no_pfam", "w") as output:
            for cluster in results:
                if cluster != None:
                    output.write(cluster)

