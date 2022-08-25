#!/usr/bin/env python3

import sys
import pandas as pd
import os
import logging
import numpy as np
from Bio import SeqIO, Seq
import argparse
import time
import csv

def mpWindowCount2(entries):
    return entries[4] + [entries[2] + 1, entries[3]] + [entries[0].count_overlap(a) for a in entries[1]]

parser = argparse.ArgumentParser(description='Spectra genetic profiling. Counting, processing, and visualization of 3-mers')
parser.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parser.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input file type', default='fasta')
parser.add_argument('-w', '--width', dest='width', type=int, help='Window width', default='3000')
parser.add_argument('-s', '--spacing', dest='spacing', type=int, help='Window spacing', default='3000')
parser.add_argument('-o', '--output', dest='output', type=str, help='Output tsv file', default='spectra_report.tsv')
parser.add_argument('-l', '--libraries', dest='libraries', action='store_true', help='Sequence names include multiple libraries, prefixed by LIBRARY_', default=False)
parser.add_argument('-p', '--proportions', dest='proportions', action='store_true', help='Return Spectra 3-mer proportions instead of raw counts', default=False)
args = parser.parse_args()

startTime = time.time()
if not os.path.exists(args.input_sequence):
    logging.error(f"Couldn't find input file '{args.input_sequence}'")
    exit()

sequences = {}
try:
    sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))
    if len(sequences.keys()) == 0:
        logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
        exit()
except ValueError:
    logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
    exit()

with open(args.output, 'w', newline='') as fileOutput:
    tsvWriter = csv.writer(fileOutput, delimiter='\t')

    queries = ["A", "C", "G", "T"]

    tsvWriter.writerow(["Library", "Sequence", "Start", "End"] + queries)

    for sequence in sequences:
        headers = sequence.split("_") if args.libraries else [os.path.basename(args.input_sequence), sequence]
        toProcess = [[sequences[sequence][i:i+args.width].seq.upper(), queries, i, i + args.width, headers] for i in range(0, len(sequences[sequence])+args.spacing, args.spacing)]
        rows = map(mpWindowCount2, toProcess)
        #rows = [mpWindowCount(sequences[sequence].seq.upper()[i:i + args.width], queries, i, i + args.width, headers) for i in range(0, len(sequences[sequence]), args.spacing)]
        tsvWriter.writerows(rows)
#            for i in range(0, len(sequences[sequence]), args.spacing):
#                rowText = mpWindowCount(sequences[sequence].seq.upper()[i:i + args.width], queries, i, i + args.width, headers)
#                tsvWriter.writerow(rowText)

        logging.info(f"Sequence {sequence} windows written to output file")
logging.info(f'Execution time in seconds: {time.time() - startTime}')