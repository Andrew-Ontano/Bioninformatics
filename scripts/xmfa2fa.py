#!/usr/bin/env python3

from Bio import AlignIO
import gfapy
import sys
import argparse

parser = argparse.ArgumentParser(description='Import mauve xmfa and output each alignment to separate file')
parser.add_argument('-i', dest='input_xmfa', type=str, required=True, help='Input xmfa file')
parser.add_argument('-o', dest='output_fasta', type=str, help='Output file prefix', default='mauve_alignments.fa')
args = parser.parse_args()

#Parse inputs as list of MSAs
inputAlign = AlignIO.parse(args.input_xmfa,"mauve")
alignments = list(inputAlign)

with open(args.output_fasta, "w") as f:
    for aln in alignments:
        for record in aln:
            f.write(">" + record.id.split('\\')[-1].replace(" ","_") + "\n")
            f.write(str(record.seq) + "\n")