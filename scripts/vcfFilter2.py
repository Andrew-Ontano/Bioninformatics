#!/usr/bin/env python3

# given an input vcf file, output a subset based on genotype filtering

import logging
import argparse
import os
import vcf
import itertools
import pandas as pd
import numpy as np
import math
import re
import gzip

# Set up logger
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='vcfFilter: a tool for processing minigraph-cactus vcf results by genotype information')
parser.add_argument('-i', '--input', dest='input_vcf', type=str, help='Input VCF file', required=True)
parser.add_argument('-o', '--output', dest='output_vcf', type=str, help='Output VCF file', required=True)
parser.add_argument('-s', '--size', dest='size', type=int, help='Number of GT columns', required=True)
parser.add_argument('-m', '--missing-allowed', dest='missing', type=int, help="Maximum number of missing alleles per row", default=3)
parser.add_argument('-g', '--groups', dest='groups', type=str, help="Sample groups to treat as same samples. Format: (A,B),(C,D),(E)", default=None)
#parser.add_argument('-t', '--targets', dest='target_genotypes', type=str, help='Genotype indices, separated by comma', nargs='*', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_vcf):
    logger.error(f"Couldn't find input VCF '{args.input_vcf}'")
    exit()

###### Config section: variables for tinkering if tinkering is needed
window_size = 10000
##### End config


if args.input_vcf.endswith('.gz'):
    with gzip.open(args.input_vcf, 'rt') as inFile:
        vcfReport = vcf.Reader(inFile)
else:
    vcfReport = vcf.Reader(open(args.input_vcf), 'r')

# Filter out results with non-calls
# Filter out only events with length change

# Build the vcf table structure
sampleNames = vcfReport.samples
if args.missing > len(sampleNames):
    sample_floor = len(sampleNames) # Minimum number of alleles to accept a vcf record with missing data
else:
    sample_floor = args.missing

outputDict = {a:{b:{c:0 for c in sampleNames} for b in list(range(0,vcfReport.contigs[a].length,window_size))} for a in vcfReport.contigs}
baselineDict = {a:{b:0 for b in list(range(0,vcfReport.contigs[a].length,window_size))} for a in vcfReport.contigs}

if args.groups is not None:
    groups = [group.split(',') for group in re.findall(r'\((.*?)\)', args.groups)]
    for record in vcfReport:
        if record.num_called >= sample_floor:
            sampleCallDict = {a.sample: a.gt_bases for a in record.samples if a.gt_bases is not None}
            sampleCallList = [a.gt_bases for a in record.samples if a.gt_bases is not None]
            window = int(math.floor(record.POS / window_size)) * window_size
            baselineDict[record.CHROM][window] += 1

            outString = record.CHROM + "\t" + str(window) + "\t" + '\t'.join([','.join(list(set(sampleCallDict[a] for a in b))) for b in groups])
            print(outString)

else:
    for record in vcfReport:
        if record.num_called >= sample_floor:
            sampleCallDict = {a.sample: a.gt_bases for a in record.samples if a.gt_bases is not None}
            sampleCallList = [a.gt_bases for a in record.samples if a.gt_bases is not None]
            window = int(math.floor(record.POS / window_size)) * window_size
            baselineDict[record.CHROM][window] += 1
            for sample in sampleCallDict:
                if len([a for a in sampleCallList if a == sampleCallDict[sample]]) == 1:
                    outputDict[record.CHROM][window][sample] += 1


    with open(f"{args.output_vcf}.tsv", 'w') as output:
        header = '\t'.join([a for a in sampleNames])
        output.write(f"Chrom\tWindow\t{header}\n")
        for chrom in outputDict.keys():
            for window in outputDict[chrom].keys():
                tallies = "\t".join([str(outputDict[chrom][window][a]) for a in outputDict[chrom][window].keys()])
                output.write(f"{chrom}\t{window}\t{tallies}\n")
    with open(f"{args.output_vcf}.index", "w") as output:
        output.write(f"Chrom\tWindow\tBaseline\n")
        for chrom in baselineDict.keys():
            for window in baselineDict[chrom].keys():
                output.write(f"{chrom}\t{window}\t{str(baselineDict[chrom][window])}\n")


#    columns = list(range(args.size))
#    targetSets = [[int(j) for j in i.split(',')] for i in args.target_genotypes]
#    targetKey = [i[0] for i in targetSets]
#    targets = list(itertools.chain.from_iterable(targetSets))
#
#    output.write(f"Chromosome\tPosition\tHit\n")
#    for record in vcfReport:
#        gts = [record.samples[i].gt_bases for i in targets]
#        gtsSets = [[record.samples[j].gt_bases for j in i] for i in targetSets]
#        # make sure no Ns present
#        if None in gts:
#            pass
#        elif "N" not in "".join(gts) and len("".join(gts)) == len(targets):
#            gtsKey = [record.samples[i].gt_bases for i in targetKey]
#            # check if both alleles the same or neither matches reference
#            if len(set(gts)) > 1:
#            #if sum([len(set(a))-1 for a in gtsSets]) == 0 and len(set(gtsKey)) == len(targetSets):
#                output.write(f"{record.CHROM}\t{record.POS}\t1\n")
#            else:
#                output.write(f"{record.CHROM}\t{record.POS}\t0\n")

