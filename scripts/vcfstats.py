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
parser.add_argument('-r', '--reference', dest='reference', type=int, help='Reference column', required=True)
parser.add_argument('-s', '--size', dest='size', type=int, help='Number of GT columns', required=True)
#parser.add_argument('-t', '--targets', dest='target_genotypes', type=str, help='Genotype names, separated by commas', required=True)
#parser.add_argument('-x', '--outgroup', dest='outgroup', type=str, help='Outgroup name', required=True)
parser.add_argument('-m', '--minimum-variant-size', dest='minimum_size', type=int, help='Minimum size of variants to include', default=1)
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
    columns.pop(args.reference)

    output.write(f"Chromosome\tPosition\tGenome\tRefSeq\tCurrSeq\tSizeDiff\n")
    for record in vcfReport:
        refSeq = record.samples[args.reference].gt_bases
        if refSeq is not None:
            if "N" not in refSeq:
                for seq in columns:
                    if record.samples[seq].gt_bases is not None:
                        if "N" not in record.samples[seq].gt_bases and abs(len(refSeq)-len(record.samples[seq].gt_bases)) >= args.minimum_size:
                            output.write(f"{record.samples[seq].site.CHROM}\t{record.samples[seq].site.POS}\t{record.samples[seq].sample}\t{refSeq}\t{record.samples[seq].gt_bases}\t{len(record.samples[seq].gt_bases)-len(refSeq)}\n")



