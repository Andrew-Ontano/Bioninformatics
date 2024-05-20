#!/usr/bin/env python3

# given an input vcf file, output a subset based on genotype filtering

import logging
import argparse
import os
import vcf
import itertools
import pandas as pd
import numpy as np

# Set up logger
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='vcfFilter: a tool for processing vcf results by genotype information')
parser.add_argument('-i', '--input', dest='input_vcf', type=str, help='Input VCF file', required=True)
parser.add_argument('-o', '--output', dest='output_vcf', type=str, help='Output VCF file', required=True)
parser.add_argument('-s', '--size', dest='size', type=int, help='Number of GT columns', required=True)
parser.add_argument('-t', '--targets', dest='target_genotypes', type=str, help='Genotype indices, separated by comma', nargs='*', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_vcf):
    logger.error(f"Couldn't find input VCF '{args.input_vcf}'")
    exit()

vcfReport = vcf.Reader(open(args.input_vcf), 'r')

# Filter out results with non-calls
# Filter out only events with length change

with open(args.output_vcf, 'w') as output:
    columns = list(range(args.size))
    targetSets = [[int(j) for j in i.split(',')] for i in args.target_genotypes]
    targetKey = [i[0] for i in targetSets]
    targets = list(itertools.chain.from_iterable(targetSets))

    output.write(f"Chromosome\tPosition\tHit\n")
    for record in vcfReport:
        gts = [record.samples[i].gt_bases for i in targets]
        gtsSets = [[record.samples[j].gt_bases for j in i] for i in targetSets]
        # make sure no Ns present
        if None in gts:
            pass
        elif "N" not in "".join(gts) and len("".join(gts)) == len(targets):
            gtsKey = [record.samples[i].gt_bases for i in targetKey]
            # check if both alleles the same or neither matches reference
            if len(set(gts)) > 1:
            #if sum([len(set(a))-1 for a in gtsSets]) == 0 and len(set(gtsKey)) == len(targetSets):
                output.write(f"{record.CHROM}\t{record.POS}\t1\n")
            else:
                output.write(f"{record.CHROM}\t{record.POS}\t0\n")


