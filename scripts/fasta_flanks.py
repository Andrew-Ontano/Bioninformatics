#!/bin/python

# given a fasta file, this program finds the first N and last M positions of each sequences, and returns a fasta file with the leading and trailing sequences only
import sys
import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='file', type=str, help='Input file name',required=True)
parser.add_argument('-o', dest='outfile', type=str, help='Output file name', default=False)
parser.add_argument('-l', dest='lead', action='store_true', help='Sequence lead to retain in bp', default=1000)
parser.add_argument('-t', dest='trail', action='store_true', help='Sequence trail to retain in bp', default=1000)
args = parser.parse_args()

if not os.path.exists(args.file):
    print("File parameter is incorrect.")
    exit()

genome = SeqIO.to_dict(SeqIO.parse(args.file,'fasta'))

if args.outfile:
    with open(args.outfile,"w") as f:
        for sequence in genome:
            if args.lead > 0:
                f.write(">" + sequence + "_leading_" + str(args.lead)+"\n")
                f.write(str(genome[sequence][0:args.lead].seq)+"\n")
            if args.trail > 0:
                f.write(">" + sequence + "_trailing_" + str(args.trail)+"\n")
                f.write(str(genome[sequence][-args.trail:].seq)+"\n")
else:
    for sequence in genome:
        if args.lead > 0:
            print(">"+sequence+"_leading_"+str(args.lead))
            print(genome[sequence][0:args.lead].seq)
        if args.trail > 0:
            print(">"+sequence+"_trailing_"+str(args.trail))
            print(genome[sequence][-args.trail:].seq)