#!/usr/bin/env python3

# Set of tools for analyzing fasta raw reads

import logging
import argparse
import os
import time
from Bio import SeqIO
import pandas as pd
import re

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='RawReadTools')
subparsers = parser.add_subparsers(title='RawReadTools', dest='subparsers', required=True)

parserLengths = subparsers.add_parser("length", description="Get list of raw read lengths")
parserLengths.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parserLengths.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input sequence format', default='fasta')
parserLengths.add_argument('-m', '--memory', dest='memory', action='store_true', help='Use memory-conservation mode', default=False)
parserLengths.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=True)

parserMask = subparsers.add_parser("missing", description="Tabulate missing data as a proportion of the length of the read")
parserMask.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parserMask.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input sequence format', default='fasta')
parserMask.add_argument('-o', '--output', dest='output_table', type=str, help='Output table file', required=True)
parserMask.add_argument('-m', '--memory', dest='memory', action='store_true', help='Use memory-conservation mode', default=False)
parserMask.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)


parserHead = subparsers.add_parser("head", description="Generate the top-N reads with missing data by proportion and by raw missing data")
parserHead.add_argument('-i', '--input', dest='input_table', type=str, help='Input tsv file', required=True)
parserHead.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

args = parser.parse_args()
startTime = time.time()

if args.verbose:
    logger.setLevel(logging.INFO)

if args.subparsers == 'missing':
    if not os.path.exists(args.input_sequence):
        logging.error(f"Couldn't find input file '{args.input_sequence}'")
        exit()
    sequences = {}
    if args.memory:
        try:
            sequences = SeqIO.index(args.input_sequence, args.sequence_format)
            if len(sequences.keys()) == 0:
                logging.error(
                    f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
                exit()
        except ValueError:
            logging.error(
                f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
            exit()
    else:
        try:
            sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))
            if len(sequences.keys()) == 0:
                logging.error(
                    f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
                exit()
        except ValueError:
            logging.error(
                f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
            exit()

    with open(args.output_table, 'w') as fileOutput:
        fileOutput.write(f'Read\tLength\tN-proportion\n')
        for seq in sequences:
            fileOutput.write(f'{seq}\t{len(sequences[seq])}\t{sequences[seq].seq.count("N")}\n')
elif args.subparsers == 'head':
    if not os.path.exists(args.input_table):
        logging.error(f"Couldn't find input file '{args.input_table}'")
        exit()
    pass
    mostMissing = []
    highProportion = []
    dataframe = pd.read_csv(args.input_table, sep='\t')
    dataframe['Proportion'] = dataframe['N-proportion']/dataframe['Length']
    subset = dataframe.loc[dataframe['Proportion'] >= 1].loc[dataframe['N-proportion'] >= 10000]
    for a in subset.iterrows():
        print(a[1]['Read'])
elif args.subparsers == 'length':
    if not os.path.exists(args.input_sequence):
        logging.error(f"Couldn't find input file '{args.input_sequence}'")
        exit()
    sequences = {}
    if args.memory:
        try:
            sequences = SeqIO.index(args.input_sequence, args.sequence_format)
            if len(sequences.keys()) == 0:
                logging.error(
                    f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
                exit()
        except ValueError:
            logging.error(
                f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
            exit()
    else:
        try:
            sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))
            if len(sequences.keys()) == 0:
                logging.error(
                    f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
                exit()
        except ValueError:
            logging.error(
                f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
            exit()
    kmers = [0] * 20
    for sequence in sequences:
        seqLen = len(sequences[sequence])
        kmers[0] += seqLen
        kmers[1] += seqLen - 1
        kmers[2] += seqLen - 2
        kmers[3] += seqLen - 3
        kmers[4] += seqLen - 4
        kmers[5] += seqLen - 5
        kmers[6] += seqLen - 6
        kmers[7] += seqLen - 7
        kmers[8] += seqLen - 8
        kmers[9] += seqLen - 9
        kmers[10] += seqLen - 10
        kmers[11] += seqLen - 11
        kmers[12] += seqLen - 12
        kmers[13] += seqLen - 13
        kmers[14] += seqLen - 14
        kmers[15] += seqLen - 15
        kmers[16] += seqLen - 16
        kmers[17] += seqLen - 17
        kmers[18] += seqLen - 18
        kmers[19] += seqLen - 19
    print(kmers)
