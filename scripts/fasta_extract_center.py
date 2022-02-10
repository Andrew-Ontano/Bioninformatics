#!/usr/bin/python3
import sys
import os
from random import shuffle
import argparse
from Bio import SeqIO
import csv
from math import ceil, floor

parser = argparse.ArgumentParser(description='Count the prevalence of a spcific k-mer across windows of a fasta file')
parser.add_argument('-i', dest='genome_fasta', type=str, help='Input genome fasta',required=True)
parser.add_argument('-l', dest='minimum_length', type=int, help='Minimum size of sequence to extract from', default='5000')
parser.add_argument('-s', dest='fragment_size', type=int, help='Size of fragment to keep', default='5000')
parser.add_argument('-o', dest='output_file', type=str, help='Output fasta file',required=False)
args = parser.parse_args()

if not os.path.exists(args.genome_fasta):
    print("Insufficient fasta file")
    exit()

genome = SeqIO.to_dict(SeqIO.parse(args.genome_fasta,'fasta'))

if args.output_file:
    outFile = open(args.output_file,'w')
    outWriter = SeqIO.FastaIO.FastaWriter(outFile)
for entry in genome:
    length = len(genome[entry])
    if length >= args.minimum_length:
        sequence = genome[entry][floor(length/2)-floor(args.fragment_size/2):floor(length/2)+ceil(args.fragment_size/2)]
        if args.output_file:
            outWriter.write_record(sequence)
        else:
            print(">"+sequence.name)
            print(sequence.seq)
if args.output_file:
    outFile.close()