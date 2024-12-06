#!/usr/bin/env python3

# given an input vcf file, output a subset based on genotype filtering

import logging
import argparse
import os

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import gzip
import re

# Set up logger
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

def readVCF(vcfFile):
    if vcfFile.endswith('.gz'):
        with gzip.open(vcfFile, 'rt') as inFile:
            for line in inFile:
                if line.startswith('#CHROM'):
                    header = [x for x in line.split('\t')]
                    break
    else:
        with open(vcfFile, 'r') as inFile:
            for line in inFile:
                if line.startswith('#CHROM'):
                    header = [x.strip() for x in line.split('\t')]
                    break

    return pd.read_csv(vcfFile, dtype=str, sep='\t', header=None, names=header, comment="#")

def filterVCF(tempVcfReport, missingAllowed):
    gtColumns = tempVcfReport.columns[9:]
    missingCounts = (tempVcfReport[gtColumns] == ".").sum(axis=1)
    return tempVcfReport[missingCounts <= missingAllowed]

def windowScheme(length, size, spacing):
    currentPosition = 0
    outputIndices = []
    while currentPosition <= length:
        outputIndices.append((currentPosition, currentPosition + size if currentPosition + size < length else length))
        currentPosition += spacing + size
        if currentPosition >= length:
            break
    return outputIndices

def pdToAlignment(df, ref=False):
    columnNames = df.columns[9:]
    nucleotideDict = {a: "" for a in columnNames}
    if ref:
        refSeq = ""
        for row, record in df.iterrows():
            alleles = [record['REF']] + record['ALT'].split(',')
            for column in columnNames:
                allele = alleles[int(record[column])] if record[column] != '.' else '-'
                nucleotideDict[column] += allele
            refSeq += alleles[0]
        referenceSeqRecord = SeqRecord(Seq(refSeq), id="Reference", name="Reference")
        alignment = MultipleSeqAlignment(records=[referenceSeqRecord])
        for column in columnNames:
            columnSeqRecord = SeqRecord(Seq(nucleotideDict[column]), id=column, name=column)
            alignment.append(columnSeqRecord)
    else:
        for row, record in df.iterrows():
            alleles = [record['REF']] + record['ALT'].split(',')
            for column in columnNames:
                allele = alleles[int(record[column])] if record[column] != '.' else '-'
                nucleotideDict[column] += allele
        alignment = MultipleSeqAlignment(records=[])
        for column in columnNames:
            columnSeqRecord = SeqRecord(Seq(nucleotideDict[column]), id=column, name=column)
            alignment.append(columnSeqRecord)
    return alignment

def buildTree(alignment):
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'nj')
    return constructor.build_tree(alignment)

def beautifyTree(tree):
    return re.sub(r'Inner\d+', '', str(tree.format(fmt='newick')).strip())

# Set up top level module argparser
parser = argparse.ArgumentParser(description='vcf2al: a tool for producing fasta alignments from vcf of SNPs')
parser.add_argument('-i', '--input', dest='input_vcf', type=str, help='Input VCF file', required=True)
#parser.add_argument('-o', '--output', dest='output_prefix', type=str, help='Output alignment prefix', default=None)
parser.add_argument('-w', '--window-size', dest='size', type=int, help='Window size in # of VCF records to generate alignments from', default=10000)
parser.add_argument('-s', '--window-overlap', dest='spacing', type=int, help='Window overlap in # of VCF records between each window. Can be negative for spaced windows', default=0)
parser.add_argument('-m', '--missing-allowed', dest='missing', type=int, help="Maximum number of missing alleles per row", default=0)
#parser.add_argument('-c', '--compare', dest='compare', action='store_true', help='Perform comparison between species with windowed trees, then report distance from species tree', default=False)
parser.add_argument('-r', '--reference', dest='reference', action='store_true', help='Incorporate reference lineage into alignments', default=False)
#parser.add_argument('-t', '--tree-file', dest='tree', type=str, help="Species tree in newick format. If absent, species tree is generated from VCF", default=None)
#parser.add_argument('-t', '--targets', dest='target_genotypes', type=str, help='Genotype indices, separated by comma', nargs='*', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-d', '--distance', dest='distance', type=str, help='Comma-separated list for genotypes for tip-distances. Calculate permutations between these samples. Use "All" to calculate all possible permutations.', default=None)
parser.add_argument('-a', '--aberrant_threshold', dest='aberrant', type=float, help='Weighted branch length threshold to call a branch aberrant. If not in distance mode, prints windows with aberrant branches', default=None)
args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_vcf):
    logger.error(f"Couldn't find input VCF '{args.input_vcf}. Exiting.'")
    exit()

# Process vcf then filter to the missing data threshold
vcfReport = readVCF(args.input_vcf)
filteredVCF = filterVCF(vcfReport, args.missing)

# Construct the full alignment form allele data

if args.compare:
# Construct a neighbor-joining tree from alignment or supplied from file
    fullBranchLength = 0
    if args.tree:
        if not os.path.exists(args.input_vcf):
            logger.error(f"Couldn't find input tree '{args.tree}. Generating a tree instead.'")
            fullAlignment = pdToAlignment(filteredVCF)
            fullTree = buildTree(fullAlignment)
            fullBranchLength = fullTree.total_branch_length()
        else:
        # get the tree from file
            fullTree = Phylo.read(args.tree, "newick")
            fullBranchLength = fullTree.total_branch_length()
    else:
        fullAlignment = pdToAlignment(filteredVCF)
        if len(fullAlignment[0]) > 0:
            # Construct NJ tree
            fullTree = buildTree(fullAlignment)
            fullBranchLength = fullTree.total_branch_length()
        else:
            logger.error(f"Insufficient variant information after filtering. Check the input VCF: {args.input_vcf}. Exiting.")
            exit()
# Prepare output headers

# Process vcf and window across chromosomes, finding alignments that match the thresholds
if args.distance:
    if args.reference:
        fullNames = list(filteredVCF.columns[8:])
        fullNames[0] = 'Reference'
    else:
        fullNames = list(filteredVCF.columns[9:])
    sampleNames = fullNames if args.distance == "All" else [a.strip() for a in args.distance.split(',') if a in fullNames]
    pairingList = [(a.strip(), b.strip()) for index, a in enumerate(sampleNames) for b in sampleNames[index + 1:]]
else:
    pairingList = []
    sampleNames = []

# Write and process header
currentRow = f"Chromosome\tStart\tEnd\tBP\tNewick\t"+"\t"
if args.distance:
    currentRow += "\t".join([f"{a[0]}-{a[1]}" for a in pairingList])
print(currentRow)

for name, group in filteredVCF.groupby('#CHROM'):
    schemes = windowScheme(len(group), args.size, args.spacing)
    group = group.reset_index().drop(columns='index')
    for scheme in schemes:
        schemeAlign = pdToAlignment(group[scheme[0]:scheme[1]], args.reference)
        schemeTree = buildTree(schemeAlign)
        schemeBranchLength = schemeTree.total_branch_length()

        if args.aberrant:
            if len([a.branch_length / schemeBranchLength for a in schemeTree.depths() if
                    a.branch_length / schemeBranchLength >= args.aberrant]) > 0:
                pass

        currentRow = f"{name}\t{group.POS[int(scheme[0])]}\t{group.POS[int(scheme[1]) - 1]}\t{sum([len(a.seq.replace('-', '')) for a in schemeAlign])}\t{beautifyTree(schemeTree)}\t"
        if args.distance:
            currentRow += "\t".join([str(schemeTree.distance(a[0], a[1])/schemeBranchLength) for a in pairingList])
        print(currentRow)
