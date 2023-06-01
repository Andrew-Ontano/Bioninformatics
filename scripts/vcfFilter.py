#!/usr/bin/env python3

# given an input vcf file, output a subset based on genotype filtering

import logging
import argparse
import os

# Set up logger
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='vcfFilter: a tool for processing vcf results by genotype information')
parser.add_argument('-i', '--input', dest='input_vcf', type=str, help='Input VCF file', required=True)
parser.add_argument('-o', '--output', dest='output_vcf', type=str, help='Output VCF file', required=True)
parser.add_argument('-t', '--targets', dest='target_genotypes', type=str, help='Genotype names, separated by commas', required=True)
parser.add_argument('-m', '--minimum-variant-size', dest='minimum_size', type=int, help='Minimum size of variants to include', default=1)
parser.add_argument('-x', '--outgroup', dest='outgroup', type=str, help='Outgroup name', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_vcf):
    logging.error(f"Couldn't find input VCF '{args.input_vcf}'")
    exit()

# Process the desired genotypes
genotypeNames = {}
outgroupNames = (None, None)

with open(args.input_vcf, 'r') as vcfInput, open(args.output_vcf, 'w') as vcfOutput:
    logger.info(f"VCF '{args.input_vcf} opened. Working through file.'")
    for line in vcfInput:
        if genotypeNames:
            lineSplit = line.strip().split('\t')
            # Process the vcf data rows if the header row has already been passed and processed
            # Set the genotype variants for the current loci
            genotypeVariants = [lineSplit[genotypeNames[a]] if genotypeNames[a] != 0 else '0' for a in genotypeNames]
            outgroupVariants = lineSplit[outgroupNames[1] if outgroupNames != 0 else '0']

            # Check if the relevant columns:
            # A) have variants that were called
            # B) are not all the same variant
            # C) are polarizable against the outgroup
            if '.' in genotypeVariants or outgroupVariants == '.' or len(set(genotypeVariants)) == 1 or outgroupVariants not in genotypeVariants:
                pass
            else:
                variants = [lineSplit[3]] + lineSplit[4].split(',')
                if sum([1 if len(variants[int(a)]) >= args.minimum_size else 0 for a in list(set(genotypeVariants + [outgroupVariants]))]) > 0:
                    vcfOutput.write(line)
        else:
            vcfOutput.write(line)
            if line.startswith('#CHROM'):
                lineSplit = line.strip().split('\t')
                for name in args.target_genotypes.strip().split(','):
                    if name in lineSplit:
                        genotypeNames[name] = lineSplit.index(name)
                    else:
                        genotypeNames[name] = 0
                if args.outgroup in lineSplit:
                    outgroupNames = (args.outgroup, lineSplit.index(args.outgroup))
                else:
                    outgroupNames = (0, args.outgroup)
                if args.outgroup in genotypeNames:
                    genotypeNames.pop(args.outgroup)
                logger.info(f"Genotype names and info processed.\nGenotypes (Dictionary of name / column indices):{genotypeNames}; Outgroup (tuple of name / column index): {outgroupNames}.\nBeginning variant processing...")
    logger.info(f"VCF '{args.input_vcf}' has been fully processed. Filtered VCF written to '{args.output_vcf}'")
