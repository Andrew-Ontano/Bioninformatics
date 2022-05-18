#!/bin/python3
import os
import argparse
from Bio import SeqIO
import csv
from math import ceil, floor

parser = argparse.ArgumentParser(description='Analyze a spectra 3-mer output for congruent areas')
parser.add_argument('-i', dest='input_gff', type=str, help='Input gff file', required=True)
parser.add_argument('-o', dest='tsv_output', type=str, help='Output tsv name', default='gff_tsv.tsv')
parser.add_argument('-q', dest='quiet', action='store_false', help='Print info on the run', default=True)
parser.add_argument('-n', dest='names', type=str, help='Criteria to save, separated by commas', default=None)
args = parser.parse_args()

if not os.path.exists(args.input_gff):
    print("Error: input file not found")
    exit()

with open(args.input_gff, 'r') as f, open(args.tsv_output, 'w') as w:
    if args.names:
        names = args.names.split(',')
        for line in f:
            if line[0] != "#":
                record = line.split("\t")
                if record[2] in names:
                    w.write(line)
    else:
        for line in f:
            if line[0] != "#":
                record = line.split("\t")
                w.write(line)
