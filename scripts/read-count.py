#!/bin/python
import sys
import os
from random import shuffle
import argparse
from Bio import SeqIO
import csv
from math import ceil, floor

parser = argparse.ArgumentParser(description='Produce table of raw read lengths')
parser.add_argument('-i', dest='genome_fastq', type=str, help='Input genome fastq')
parser.add_argument('-o', dest='report_output', type=str, help='Output file name', default='repeat_report.tsv')
args = parser.parse_args()
if not args.genome_fastq:
    print("Please specify the proper input file")
    exit()

# make sure both file paths exist
if not os.path.exists(args.genome_fastq):
    print("File parameter is incorrect. Usage: 'python read-counts.py -i $genome_fasta'")
    exit()

# load genome file
genome = SeqIO.to_dict(SeqIO.parse(args.genome_fastq,'fastq'))
# loop through remaining entries in order to find

readLengths={}

for rawRead in genome:
    rawReadLength = len(genome[rawRead].seq)
    if rawReadLength in readLengths:
        readLengths[rawReadLength] += 1
    else:
        readLengths[rawReadLength] = 1

output = open(args.report_output, 'w', newline='')
tsvWriter = csv.writer(output, delimiter='\t')
tsvWriter.writerow(["Length", "Reads"])

for readLength in readLengths:
    rowText=[readLength,readLengths[readLength]]
    tsvWriter.writerow(rowText)
output.close()
