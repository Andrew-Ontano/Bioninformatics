#!/usr/bin/env python3

# given an input vcf file, output a subset based on genotype filtering

import logging
import argparse
import os
import pandas as pd
import numpy as np
from math import floor, ceil
# Set up logger
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='vcfFilter: a tool for processing vcf results by genotype information')
parser.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input ts file', required=True)
parser.add_argument('-o', '--output', dest='output_tsv', type=str, help='Output tsv file', required=True)
parser.add_argument('-s', '--size', dest='size', type=int, help='Bin size', default=2000000)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_tsv):
    logger.error(f"Couldn't find input TSV '{args.input_tsv}'")
    exit()

df = pd.read_csv(args.input_tsv, sep='\t')
df_grouped = df.groupby("Chromosome")

with open(args.output_tsv, 'w') as output:
    output.write(f"Chromosome\tStart\tEnd\tHetSites\tTotalSites\n")
    for name, group in df_grouped:
        bins = ceil(max(group["Position"])/args.size)
        for groupbin in range(bins):
            width = groupbin * args.size
            groupBinned = group[group['Position'].between(width+1, width+args.size)]
            output.write(f"{name}\t{width+1}\t{width+args.size}\t{sum(groupBinned['Hit'])}\t{groupBinned.shape[0]}\n")
