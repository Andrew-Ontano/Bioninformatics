#!/usr/bin/env python3

# Given a fasta alignment file and a gapless Spectra, output gapless fasta file

import sys
import pandas as pd
import os
import logging
import numpy as np
from Bio import SeqIO, Seq

logging.basicConfig(level=logging.INFO)

if len(sys.argv) < 3:
    logging.warning('Insufficient files supplied')
    exit()
if not os.path.exists(sys.argv[1]) or not os.path.exists(sys.argv[2]):
    logging.warning('One or both files do not exist.')
    exit()

sequences = {}
try:
    sequences = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))
    if len(sequences.keys()) == 0:
        logging.error(
            f"Sequence file '{sys.argv[1]}' could not be loaded in fasta format or has incorrectly formatted sequences."
        )
        exit()
except ValueError:
    logging.error(f"Sequence file '{sys.argv[1]}' could not be loaded in fasta format.")
    exit()

with open(sys.argv[2]) as inFile:
    data = pd.read_csv(inFile, sep='\t')

outSequences = {seq: Seq.Seq('') for seq in sequences}

for start, end in zip(data['Start'].unique(), data['End'].unique()):
    for seq in sequences:
        outSequences[seq] += sequences[seq][start:end]

for outSequence in outSequences:
    outSequences[outSequence].seq = outSequences[outSequence].seq.ungap()
    with open(f'{outSequence}.fna', 'w') as outFile:
        SeqIO.write(outSequences[outSequence], outFile, 'fasta')

