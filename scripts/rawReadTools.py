#!/usr/bin/env python3

# Set of tools for analyzing fasta raw reads

import logging
import argparse
import os
import time
from Bio import SeqIO, Seq
import re

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='RawReadTools')
subparsers = parser.add_subparsers(title='RawReadTools', dest='subparsers', required=True)

parserMask = subparsers.add_parser("missing", description="Tabulate missing data as a proportion of the length of the read")
parserMask.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parserMask.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input sequence format', default='fasta')
parserMask.add_argument('-o', '--output', dest='output_table', type=str, help='Output table file', required=True)
parserMask.add_argument('-m', '--memory', dest='memory', action='store_true', help='Use memory-conservation mode', default=False)
parserMask.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

startTime = time.time()
if not os.path.exists(args.input_sequence):
    logging.error(f"Couldn't find input file '{args.input_sequence}'")
    exit()


sequences = {}
if args.memory:
    try:
        sequences = SeqIO.index(args.input_sequence, args.sequence_format)
        if len(sequences.keys()) == 0:
            logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
            exit()
    except ValueError:
        logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
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

#SeqIO.write((sequences[a] for a in sequences), args.output_table, args.sequence_format)
