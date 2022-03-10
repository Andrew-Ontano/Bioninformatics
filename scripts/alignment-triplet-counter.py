#!/bin/python
import sys
import os
from random import shuffle
import argparse
from Bio import SeqIO
import csv
from math import ceil, floor

parser = argparse.ArgumentParser(description='Count the prevalence of 3-mers from an alignment')
parser.add_argument('-i', dest='genome_fasta', type=str, help='Input genome fasta')
parser.add_argument('-w', dest='window_width', type=int, help='Window width', default='1000')
parser.add_argument('-s', dest='window_spacing', type=int, help='Window spacing', default='1000')
parser.add_argument('-o', dest='report_output', type=str, help='Output file name', default='repeat_report.tsv')
args = parser.parse_args()
if not args.genome_fasta:
    print("Please specify the proper input file")
    exit()

# make sure both file paths exist
if not os.path.exists(args.genome_fasta):
    print("File parameter is incorrect. Usage: 'python repeat-windows.py -i $genome_fasta -q '$query_sequence'")
    exit()

# load genome file
genome = SeqIO.to_dict(SeqIO.parse(args.genome_fasta,'fasta'))
# loop through remaining entries in order to find
output = open(args.report_output, 'w',newline='')
tsvWriter = csv.writer(output, delimiter='\t')

window_lead = floor(args.window_width/2)
window_tail = ceil(args.window_width/2)


#bases = ["A","C","G","T","N","-"]
bases = ["A","C","G","T"]
queries = [bases[a]+bases[b]+bases[c] for c in range(0,len(bases)) for b in range(0,len(bases)) for a in range(0,len(bases))]

tsvWriter.writerow(["Scaffold", "Start", "End"]+queries)

for entry in genome:
    i = window_lead
    while i+window_tail < len(genome[entry]) :
        rowSeq = str(genome[entry][i-window_lead:i+window_tail].seq).upper()
        rowCounts = [str(rowSeq.count(query)) for query in queries]
        rowText = [
            entry,
            i-window_lead+1,
            i+window_tail
        ]
        rowText += rowCounts
        tsvWriter.writerow(rowText)
        i += args.window_spacing
    rowSeq = str(genome[entry][i-window_lead:i+window_tail].seq).upper()
    rowCounts = [str(rowSeq.count(query)) for query in queries]
    rowText = [
            entry,
            i-window_lead+1,
            len(genome[entry])
    ]
    rowText += rowCounts
    tsvWriter.writerow(rowText)
output.close()
