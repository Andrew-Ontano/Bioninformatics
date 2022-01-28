#!/bin/python3
import sys
import os
from random import shuffle
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Extend the flanking region around BUSCOs outward in both directions')
parser.add_argument('-g', dest='genome_fasta', type=str, help='Input genome fasta')
parser.add_argument('-b', dest='busco_folder', type=str, help='Input busco run folder')
parser.add_argument('-r', dest='random_subsamples', type=int, help='Random subsamples to take', default='0')
parser.add_argument('-f', dest='flanking_length', type=int, help='Flanking region in both directions', default='100000')
parser.add_argument('-o', dest='fasta_output', type=str, help='Output file name', default='extended_buscos.fasta')
args = parser.parse_args()
if not args.genome_fasta or not args.busco_folder:
    print("Please specify the proper input files")
    exit()

# make sure both file paths exist
if not os.path.exists(args.genome_fasta) or not os.path.exists(args.busco_folder):
    print("One or both file parameters are incorrect. Usage: 'python extrapolate_regions.py $genome_fasta $busco_result_path [$random_subsamples]'")
    exit()

tablePath = args.busco_folder+"/"+[f for f in os.listdir(args.busco_folder) if f.startswith('full_table')][0]

# read table of buscos, and remove missing entries as well as headers
with open(tablePath) as f:
    tableContents = [a.strip() for a in f.readlines() if "Missing" not in a and "#" not in a]

# prune list to just N random entries if $random_subsamples supplied
if args.random_subsamples:
    listSize = len(tableContents)
    if args.random_subsamples < listSize:
        shuffle(tableContents)
        tableContents = tableContents[0:args.random_subsamples]

# load genome file
genome = SeqIO.to_dict(SeqIO.parse(args.genome_fasta,'fasta'))
# loop through remaining entries in order to find
fastaOutput = open(args.fasta_output, 'w')
fastaWriter = SeqIO.FastaIO.FastaWriter(fastaOutput)
for entry in tableContents:
    values = entry.split('\t')
    contigSize = len(genome[values[2]])
    left_flank = (int(values[3]) - args.flanking_length) if (int(values[3]) - args.flanking_length) >= 0 else 0
    right_flank = (int(values[4]) + args.flanking_length) if (int(values[4]) + args.flanking_length) < contigSize else contigSize

    seqRecord = genome[values[2]][left_flank:right_flank]
    seqRecord.id = values[0]+"_"+seqRecord.id
    seqRecord.description = "from position "+str(left_flank)+"-"+str(right_flank)
    fastaWriter.write_record(seqRecord)
    print(values[0], values[2], left_flank, right_flank)
fastaOutput.close()
