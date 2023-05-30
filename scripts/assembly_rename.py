#!/usr/bin/env python3
import argparse
import os
from Bio import SeqIO
import logging

# setup for cl args and logging
logging.basicConfig(level=logging.INFO)
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='file', type=str, help='Input file name.', required=True)
args = parser.parse_args()

# check for file
if not os.path.exists(args.file):
    print("File parameter is incorrect.")
    exit()

genome = SeqIO.to_dict(SeqIO.parse(args.file, 'fasta'))

print("Sequence\tLength")
for sequence in genome:
    print(sequence+"\t"+str(len(genome[sequence]))
