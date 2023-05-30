#!/usr/bin/env python3

# Tool for taking an input TRF dat file, and processing pertinent info

import logging
import argparse
import os
import time
from Bio import SeqIO, Seq
import re
import pandas as pd

# Set up logger
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='trfFilter: a tool for processing information from TRF outputs')

parser.add_argument('-i', '--input', dest='input_dat', type=str, help='Input dat file', required=True)
parser.add_argument('-o', '--output', dest='output_dat', type=str, help='Output dat file', default=None)
parser.add_argument('-l', '--limit', dest='size_limit', type=int, help='Threshold to split motif sizes', default=100)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-m', '--minimum', dest='minimum_length', type=int, help='Minimum length in bases to keep a sequence', default=1)

args = parser.parse_args()
if args.verbose:
    logger.setLevel(logging.INFO)

# Process the intitial dat file into a condensed format
outputLines=[]
currentSequence = ''
with open(args.input_dat, 'r') as fileIn:
    for line in fileIn:
        line = line.strip()
        if line.startswith("Sequence"):
            currentSequence = line.split(':')[1][1:].split(' ')[0]
        elif re.match("^[0-9]", line):
            #outputLines.append(f"{currentSequence} {line}")
            outputLines.append([currentSequence]+line.split(" ")[0:14])

# Set up pandas df
trfDF = pd.DataFrame(outputLines, columns=['sequence', 'start', 'end', 'period_size', 'number_copies', 'consensus_size', 'percent_matches', 'percent_indels', 'alignment_score', 'percent_a', 'percent_c', 'percent_g', 'percent_t', 'percent_entropy', 'consensus_sequence'])
#trfGrouped = trfDF.groupby(by='sequence')

# process each sequence to find 1) best scoring hits and 2) non-overlapping hits
#simplifiedGroups = pd.DataFrame(None, columns=['sequence', 'start', 'end', 'period_size', 'number_copies', 'consensus_size', 'percent_matches', 'percent_indels', 'alignment_score', 'percent_a', 'percent_c', 'percent_g', 'percent_t', 'percent_entropy', 'consensus_sequence', 'region_sequence'])
#for group in trfGrouped:
#    currentIndices = (0, 0)
#    # iterate through the entire group
#    for index in group[1].index:
#        # check if the length is long enough to keep
#        currentTRF = group[1].loc[index]
#        if 1 + int(currentTRF['end']) - int(currentTRF['start']) >= args.minimum_length:
#            # check if the indices are not overlapping with others
#            if int(currentTRF['start']) > currentIndices[1]:
#                currentIndices = (int(currentTRF['start']), int(currentTRF['end']))
#                simplifiedGroups.loc[len(simplifiedGroups.index)] = currentTRF

if args.output_dat:
    trfDF.to_csv(args.output_dat, sep='\t', index=False)
    #simplifiedGroups.to_csv(args.output_dat, sep='\t', index=False)
#logger.info(
#    f"Size info for {args.input_dat}- Total: {sum([len(simplifiedGroups.loc[a]['region_sequence']) for a in simplifiedGroups.index])}, <{args.size_limit}: {sum([len(simplifiedGroups.loc[a]['region_sequence']) for a in simplifiedGroups.index if len(simplifiedGroups.loc[a]['consensus_sequence']) < args.size_limit])}, >={args.size_limit}: {sum([len(simplifiedGroups.loc[a]['region_sequence']) for a in simplifiedGroups.index if len(simplifiedGroups.loc[a]['consensus_sequence']) >= args.size_limit])}"
#)
#logger.info(f"Number of enriched reads: {len(list(set([simplifiedGroups.loc[a]['sequence'] for a in simplifiedGroups.index if len(simplifiedGroups.loc[a]['region_sequence']) >= args.minimum_length])))}")
