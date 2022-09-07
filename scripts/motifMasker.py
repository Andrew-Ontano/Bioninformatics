#!/usr/bin/env python3

# Given an input motif, masks all instances of that sequence with Ns

import logging
import argparse
import os
import time
from Bio import SeqIO, Seq
import re

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

def rc(sequence):
    sequence = Seq.Seq(sequence)
    return str(sequence.reverse_complement())

# Set up top level module argparser
parser = argparse.ArgumentParser(description='MotifMasker')

parser.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parser.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input sequence format', default='fasta')
parser.add_argument('-o', '--output', dest='output_sequence', type=str, help='Output sequence file', required=True)
parser.add_argument('-q', '--query-motif', dest='motif', type=str, help='Query motif', required=True)
parser.add_argument('-m', '--memory', dest='memory', action='store_true', help='Use memory-conservation mode', default=False)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-r', '--reverse-complement', dest='reverse', action='store_true', help='Mask reverse complement of motif', default=False)

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

    removalList = f'({args.motif})|({rc(args.motif)})' if args.reverse else f'({args.motif})'
    replacementSeq = "".join("N" * len(args.motif))
    with open(args.output_sequence, 'w') as fileOutput:
        for seq in sequences:
            fileOutput.write(f'>{sequences[seq].id}\n')
            fileOutput.write(f'{re.sub(removalList, replacementSeq, str(sequences[seq].seq), flags=re.IGNORECASE)}\n')
else:
    try:
        sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))
        if len(sequences.keys()) == 0:
            logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
            exit()
    except ValueError:
        logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
        exit()

    if args.reverse:
        removalList = [args.motif, rc(args.motif)]
        replacementSeq = "".join("N" * len(args.motif))

        for seq in sequences:
            sequences[seq].seq = sequences[seq].seq.replace(removalList[0], replacementSeq).replace(removalList[1], replacementSeq)

    else:
        removalList = args.motif
        replacementSeq = "".join("N" * len(args.motif))
        for seq in sequences:
            sequences[seq].seq = sequences[seq].seq.replace(removalList, replacementSeq)
    SeqIO.write((sequences[a] for a in sequences), args.output_sequence, args.sequence_format)
