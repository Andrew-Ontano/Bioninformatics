#!/usr/bin/env python3

# given an input vcf file, output a subset based on genotype filtering

import logging
import argparse
import os
import vcf
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
parser.add_argument('-t', '--targets', dest='target_genotypes', type=str, help='Genotype indices, separated by comma', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-r', '--reference', dest='reference', type=int, help='Reference column', required=True)
parser.add_argument('-m', '--hom', dest='hom', action='store_true', help='Check for hom', default=False)

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
    targets = [int(i) for i in args.target_genotypes.split(',')]
    output.write(f"Chromosome\tPosition\tHit\n")
    for record in vcfReport:
        gts = [record.samples[i].gt_bases for i in targets]
        reference = record.samples[args.reference].gt_bases
        # make sure no Ns present
        if None in gts:
            pass
        elif "N" not in reference+"".join(gts) and len(reference) == 1 and len("".join(gts)) == len(targets):
            # check if both alleles the same or neither matches reference
            if reference in gts and not all(i == gts[0] for i in gts):
                output.write(f"{record.CHROM}\t{record.POS}\t{0 if args.hom else 1}\n")

            else:
                output.write(f"{record.CHROM}\t{record.POS}\t{1 if args.hom else 0}\n")
