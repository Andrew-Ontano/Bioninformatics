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
parser.add_argument('-g', '--genotypes', dest='genotypes_count', type=int, help='Number of genotypes', default=1)
parser.add_argument('-t', '--targets', dest='target_columns', type=str, help='Columns for comparison, separated by commas', required=True)
parser.add_argument('-m', '--minimum-variant-size', dest='minimum_size', type=int, help='Minimum size of variants to include', default=1)
parser.add_argument('-x', '--outgroup', dest='outgroup_column', type=int, help='Numeric values of outgroup. 0=reference', default=0)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_vcf):
    logging.error(f"Couldn't find input VCF '{args.input_vcf}'")
    exit()

# set a tuple with the boundaries for genotypes
genotypeColumns = [a+9 for a in range(0, args.genotypes_count)]
targetColumns = args.target_columns.split(',')
# switch to reference mode if reference included, or outgroup mode if reference is the outgroup
if '0' in targetColumns:
    genotypeTargetColumns = [genotypeColumns[int(a)-1] for a in targetColumns if a != '0']
    outgroupColumn = 9 + args.outgroup_column - 1
    with open(args.input_vcf, 'r') as vcfInput, open(args.output_vcf, 'w') as vcfOutput:
        for line in vcfInput:
            if line.startswith('#'):
                vcfOutput.write(line)
            else:
                lineSplit = line.strip().split('\t')
                genotypeValues = [lineSplit[a] for a in genotypeTargetColumns]
                outgroupValue = lineSplit[outgroupColumn]
                # make sure variants are called
                if '.' not in genotypeValues and outgroupValue != '.':
                    # add the reference ('0') to the list before making the collapsible list

                    if (outgroupValue in genotypeValues or outgroupValue == '0') and len(set(['0'] + genotypeValues)) > 1:
                        variants = [lineSplit[3]] + lineSplit[4].split(',')
                        if len([a for a in list(set(['0'] + genotypeValues)) if len(variants[int(a)]) >= args.minimum_size]) > 0:
                            #print(line)
                            vcfOutput.write(line)
elif args.outgroup_column == 0:
    logging.error("This is not supported yet, sorry!")

    #genotypeTargetColumns = [genotypeColumns[int(a)-1] for a in targetColumns if a != '0']
    #with open(args.input_vcf, 'r') as vcfInput, open(args.output_vcf, 'w') as vcfOutput:
        #for line in vcfInput:
            #if line.startswith('#'):
                #vcfOutput.write(line)
        #run comparison here
    pass
else:
    genotypeTargetColumns = [genotypeColumns[int(a)-1] for a in args.target_columns.split(',')]
    outgroupColumn = 9 + args.outgroup_column - 1
    with open(args.input_vcf, 'r') as vcfInput, open(args.output_vcf, 'w') as vcfOutput:
        for line in vcfInput:
            if line.startswith('#'):
                vcfOutput.write(line)
            else:
                lineSplit = line.strip().split('\t')
                genotypeValues = [lineSplit[a] for a in genotypeTargetColumns]
                outgroupValue = lineSplit[outgroupColumn]

                # make sure variants are called
                if '.' not in genotypeValues and outgroupValue != '.':
                    if outgroupValue in genotypeValues and len(set(genotypeValues)) > 1:
                        variants = [lineSplit[3]] + lineSplit[4].split(',')
                        if len([a for a in list(set(genotypeValues)) if len(variants[int(a)]) >= args.minimum_size]) > 0:
                            vcfOutput.write(line)
