#!/bin/python

# given a fasta file, this program finds the first N and last M positions of each sequences, and returns a fasta file with the leading and trailing sequences only
import sys
import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='file', type=str, help='Input file name',required=True)
args = parser.parse_args()

if not os.path.exists(args.file):
    print("File parameter is incorrect.")
    exit()

genome = SeqIO.to_dict(SeqIO.parse(args.file,'fasta'))

print("Sequence\tLength")
for sequence in genome:
    print(sequence+"\t"+str(len(genome[sequence])))