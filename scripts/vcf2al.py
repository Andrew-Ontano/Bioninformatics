#!/usr/bin/env python3

# given an input vcf file, output a subset based on genotype filtering

import logging
import argparse
import os

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
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
                    header = [x.strip() for x in line.split('\t')]
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

# Initiates an empty alignment with the needed names
def initAlignment(df, ref=False):
    columnNames = df.columns[9:]
    alignment = MultipleSeqAlignment(records=[])
    if ref:
        currentRecord = SeqRecord(Seq(""), id="Reference", name="Reference")
        alignment.append(currentRecord)
    for column in columnNames:
        currentRecord = SeqRecord(Seq(""), id=column, name=column)
        alignment.append(currentRecord)
    return alignment

def pdToAlignment(df, ref=False):
    alignment = initAlignment(df, ref)
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
        nucleotideDict["Reference"] = refSeq
        for a in range(len(alignment)):
            alignment[a].seq = Seq(nucleotideDict[alignment[a].name])
    else:
        for row, record in df.iterrows():
            alleles = [record['REF']] + record['ALT'].split(',')
            for column in columnNames:
                allele = alleles[int(record[column])] if record[column] != '.' else '-'
                nucleotideDict[column] += allele

        for a in range(len(alignment)):
            alignment[a].seq = Seq(nucleotideDict[alignment[a].name])

    return alignment

### Stored until we're sure this is less efficient than the redone alignments
    #     alignment = MultipleSeqAlignment(records=[referenceSeqRecord])
    #
    #     for column in columnNames:
    #         columnSeqRecord = SeqRecord(Seq(nucleotideDict[column]), id=column, name=column)
    #         alignment.append(columnSeqRecord)
    # else:
    #     for row, record in df.iterrows():
    #         alleles = [record['REF']] + record['ALT'].split(',')
    #         for column in columnNames:
    #             allele = alleles[int(record[column])] if record[column] != '.' else '-'
    #             nucleotideDict[column] += allele
    #     alignment = MultipleSeqAlignment(records=[])
    #     for column in columnNames:
    #         columnSeqRecord = SeqRecord(Seq(nucleotideDict[column]), id=column, name=column)
    #         alignment.append(columnSeqRecord)
    # return alignment

## TODO: modify method to accept multiple tree construction methods
def buildTree(alignment):
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'nj')
    return constructor.build_tree(alignment)

# Method for formatting the tree in a human-readable format (Newick)
def beautifyTree(tree):
    return re.sub(r'Inner\d+', '', str(tree.format(fmt='newick')).strip())

# Function for testing formatting issues
def testTree(tree):
    pairingList = [(a.strip(), b.strip()) for index, a in enumerate(sampleNames) for b in sampleNames[index + 1:]]


# Set up top level module argparser
parser = argparse.ArgumentParser(description='vcf2al: a tool for producing fasta alignments from vcf of SNPs')
parser.add_argument('-i', '--input', dest='input_vcf', type=str, help='Input VCF file', required=True)
parser.add_argument('-o', '--output', dest='output_tsv', type=str, help='Output alignment prefix', default='vcf2al_out.tsv')
parser.add_argument('-w', '--window-size', dest='size', type=int, help='Window size in # of VCF records to generate alignments from', default=10000)
parser.add_argument('-s', '--window-overlap', dest='spacing', type=int, help='Window overlap in # of VCF records between each window. Can be negative for spaced windows', default=0)
parser.add_argument('-m', '--missing-allowed', dest='missing', type=int, help="Maximum number of missing alleles per row", default=0)
parser.add_argument('-b', '--blacklist', dest='blacklist', type=str, help='Remove species from analysis. Separated by commas', default=None)
parser.add_argument('-c', '--compare', dest='compare', action='store_true', help='Perform comparison between species with windowed trees, then report distance from species tree', default=False)
parser.add_argument('-r', '--reference', dest='reference', action='store_true', help='Incorporate reference lineage into alignments', default=False)
parser.add_argument('-t', '--tree-file', dest='tree', type=str, help="Species tree in newick format. If absent, species tree is generated from VCF", default=None)
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

with open(args.output_tsv, 'w') as outFile:
    totalAlign = initAlignment(filteredVCF, ref=args.reference)

    currentRow = f"Type\tChromosome\tStart\tEnd\tBP\tNewick"
    if args.distance:
        currentRow += "\t" + "\t".join([f"{a[0]}-{a[1]}" for a in pairingList])
    outFile.write(currentRow+"\n")

    for name, group in filteredVCF.groupby('#CHROM'):
        schemes = windowScheme(len(group), args.size, args.spacing)
        group = group.reset_index().drop(columns='index')
        if args.compare:
            chromosomeAlign = pdToAlignment(group, args.reference)
            totalAlign += chromosomeAlign
            chromosomeTree = buildTree(chromosomeAlign)
            chromosomeBranchLength = chromosomeTree.total_branch_length()
            currentRow = f"Chromosome\t{name}\t1\t{len(chromosomeAlign[0])}\t{sum([len(a.seq.replace('-', '')) for a in chromosomeAlign])}\t{beautifyTree(chromosomeTree)}"
            if args.distance:
                currentRow += "\t" + "\t".join([str(chromosomeTree.distance(a[0], a[1])/chromosomeBranchLength) for a in pairingList])
            outFile.write(currentRow + "\n")

        for scheme in schemes:

            schemeAlign = pdToAlignment(group[scheme[0]:scheme[1]], args.reference)
            schemeTree = buildTree(schemeAlign)
            schemeBranchLength = schemeTree.total_branch_length()

            if args.aberrant:
                if len([a.branch_length / schemeBranchLength for a in schemeTree.depths() if
                        a.branch_length / schemeBranchLength >= args.aberrant]) > 0:
                    pass

            currentRow = f"Window\t{name}\t{group.POS[int(scheme[0])]}\t{group.POS[int(scheme[1]) - 1]}\t{sum([len(a.seq.replace('-', '')) for a in schemeAlign])}\t{beautifyTree(schemeTree)}"
            if args.distance:
                currentRow += "\t" + "\t".join([str(schemeTree.distance(a[0], a[1])/schemeBranchLength) for a in pairingList])
                # currentRow += "\t" + "\t".join([str(schemeTree.distance(a[0], a[1])/schemeBranchLength) for a in pairingList])
            outFile.write(currentRow+"\n")

    if args.compare:
        totalTree = buildTree(totalAlign)
        totalBranchLength = totalTree.total_branch_length()
        currentRow = f"Full\tAll\t1\t{len(totalAlign[0])}\t{sum([len(a.seq.replace('-', '')) for a in totalAlign])}\t{beautifyTree(totalTree)}"
        if args.distance:
            currentRow += "\t" + "\t".join(
                [str(totalTree.distance(a[0], a[1]) / totalBranchLength) for a in pairingList])
        outFile.write(currentRow + "\n")