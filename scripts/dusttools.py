#!/usr/bin/env python3

# Set of tools for processing DUSTmasker output fasta.

import logging
import argparse
import os
import time
from Bio import SeqIO, Seq
import re

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='DUSTtools')
subparsers = parser.add_subparsers(title='DUSTtools', dest='subparsers', required=True)

parserMask = subparsers.add_parser("mask", description="Mask either the repeat-enriched or the repeat-deficient sequence")
parserMask.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parserMask.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input sequence format', default='fasta')
parserMask.add_argument('-o', '--output', dest='output_sequence', type=str, help='Output sequence file', required=True)
parserMask.add_argument('-r', '--keep-repeats', dest='keep_repeats', action='store_true', help='Mask all sequences except repeats, instead of only repetitive elements', default=False)
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

    removalList = r'[ACGT]' if args.keep_repeats else r'[acgt]'
    with open(args.output_sequence, 'w') as fileOutput:
        for seq in sequences:
            fileOutput.write(f'>{sequences[seq].id}\n')
            fileOutput.write(f'{re.sub(removalList, "N", str(sequences[seq].seq))}\n')
else:
    try:
        sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))
        if len(sequences.keys()) == 0:
            logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
            exit()
    except ValueError:
        logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
        exit()

    removalList = ['A', 'C', 'G', 'T'] if args.keep_repeats else ['a', 'c', 'g', 't']

    for seq in sequences:
        sequences[seq].seq = sequences[seq].seq.replace(removalList[0], 'N').replace(removalList[1], 'N').replace(
            removalList[2], 'N').replace(removalList[3], 'N')
    SeqIO.write((sequences[a] for a in sequences), args.output_sequence, args.sequence_format)
