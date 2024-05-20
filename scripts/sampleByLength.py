#!/usr/bin/env python3

# Tool for taking sequences in fasta format and random sampling until they exceed a target cumulative length
import logging
import argparse
import os
import time
from Bio import SeqIO, Seq
import re
from random import shuffle

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='sampleByLength')

parser.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parser.add_argument('-o', '--output', dest='output_sequence', type=str, help='Output sequence file', required=True)
parser.add_argument('-t', '--target', dest='target_length', type=int, help='Target length in bases', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-m', '--minimum', dest='minimum_length', type=int, help='Minimum length in bases to keep a sequence', default=1)


args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, 'fasta'))
outputSequences = {}
targetLength = 0
seeds = list(range(len(sequences)))
shuffle(seeds)
while targetLength < args.target_length:
    # add next random seed
    nextSeed = seeds.pop(0)
    nextSequence = list(sequences.keys())[nextSeed]
    if len(sequences[nextSequence]) >= args.minimum_length:
        outputSequences[nextSequence] = sequences[nextSequence]
        targetLength += len(sequences[nextSequence])
SeqIO.write(outputSequences.values(), args.output_sequence, 'fasta')
