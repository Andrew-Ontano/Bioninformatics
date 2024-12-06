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

def pdToAlignment(df, ref=True):
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
    #return calculator.get_distance(alignment), constructor.build_tree(alignment)
    return constructor.build_tree(alignment)

def beautifyTree(tree):
    return re.sub(r'Inner\d+', '', str(tree.format(fmt='newick')).strip())

# Set up top level module argparser
parser = argparse.ArgumentParser(description='vcf2al: a tool for producing fasta alignments from vcf of SNPs')
parser.add_argument('-i', '--input', dest='input_vcf', type=str, help='Input VCF file', required=True)
parser.add_argument('-o', '--output', dest='output_prefix', type=str, help='Output alignment prefix', default=None)
parser.add_argument('-w', '--window-size', dest='size', type=int, help='Window size in # of VCF records to generate alignments from', default=10000)
parser.add_argument('-s', '--window-overlap', dest='spacing', type=int, help='Window overlap in # of VCF records between each window. Can be negative for spaced windows', default=0)
parser.add_argument('-m', '--missing-allowed', dest='missing', type=int, help="Maximum number of missing alleles per row", default=3)
parser.add_argument('-c', '--compare', dest='compare', action='store_true', help='Perform comparison between species with with windowed trees, then report distance from species tree', default=False)
parser.add_argument('-t', '--tree-file', dest='tree', type=str, help="Species tree in newick format. If absent, species tree is generated from VCF", default=None)
#parser.add_argument('-t', '--targets', dest='target_genotypes', type=str, help='Genotype indices, separated by comma', nargs='*', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-d', '--distance', dest='distance', type=str, help='Comma-separated list for genotypes for tip-distances. Calculate permutations between these samples. Use "All" to calculate all possible permutations.', default=None)
parser.add_argument('-a', '--aberrant_threshold', dest='aberrant', type=float, help='Weighted branch length threshold to call a branch aberrant. If not in distance mode, prints windows with aberrant branches', default=0.2)
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
if args.distance:
    fullTreeNames = [a.name for a in fullTree.get_terminals()]
    sampleNames = fullTreeNames if args.distance == "All" else [a.strip() for a in args.distance.split(',') if a in fullTreeNames]
    pairingList = [(a.strip(), b.strip()) for index, a in enumerate(sampleNames) for b in sampleNames[index + 1:]]
    print(f"Chromosome\tStart\tEnd\tBP\tNewick\t"+"\t".join([f"{a[0]}-{a[1]}" for a in pairingList]))
    print(f"All\t0\t0\t0\t{beautifyTree(fullTree)}\t"+"\t".join([str(fullTree.distance(a[0], a[1])/fullBranchLength) for a in pairingList]))
    for name, group in filteredVCF.groupby('#CHROM'):
        schemes = windowScheme(len(group), args.size, args.spacing)
        group = group.reset_index().drop(columns='index')
        for scheme in schemes:
            schemeAlign = pdToAlignment(group[scheme[0]:scheme[1]])
            schemeTree = buildTree(schemeAlign)
            schemeBranchLength = schemeTree.total_branch_length()
            print(f"{name}\t{group.POS[int(scheme[0])]}\t{group.POS[int(scheme[1]) - 1]}\t{sum([len(a.seq.replace('-', '')) for a in schemeAlign])}\t{beautifyTree(schemeTree)}\t"+"\t".join([str(schemeTree.distance(a[0], a[1])/schemeBranchLength) for a in pairingList]))
else:
    print(f"Chromosome\tStart\tEnd\tBP\tNewick")
    for name, group in filteredVCF.groupby('#CHROM'):
        schemes = windowScheme(len(group), args.size, args.spacing)
        group = group.reset_index().drop(columns='index')
        for scheme in schemes:
            schemeAlign = pdToAlignment(group[scheme[0]:scheme[1]])
            schemeTree = buildTree(schemeAlign)
            schemeBranchLength = schemeTree.total_branch_length()
            if len([a.branch_length/schemeBranchLength for a in schemeTree.depths() if a.branch_length/schemeBranchLength >= args.aberrant]) > 0:
                print(
                    f"{name}\t{group.POS[int(scheme[0])]}\t{group.POS[int(scheme[1]) - 1]}\t{sum([len(a.seq.replace('-', '')) for a in schemeAlign])}\t{beautifyTree(schemeTree)}")

# Create a nested dictionary with the structure for contigs > indices > samples ---- this is based on indices
#outputDict = {a:{b:{c:"" for c in sampleNames} for b in list(range(0,vcfReport.contigs[a].length,args.window))} for a in vcfReport.contigs}
#for record in vcfReport:
#    if record.num_called >= sample_floor:
#        window = int(math.floor(record.POS / args.window)) * args.window
#        for sample in record.samples:
#            outputDict[record.CHROM][window][sample.sample] += sample.gt_bases if sample.gt_bases is not None else "-"
#        outputDict[record.CHROM][window]["Ref"] += record.alleles[0]

# Create a nested dictionary with the structure for contigs - generated from N-bp sliding windows
# Windows = window size
#outputDict = {a:{b for b in } for a in vcfReport.contigs}

# # Print header
# outString = f"Chromosome\tIndex\tLength\tTree"
# if args.distance:
#     namesList = sampleNames
#     samplePairs = []
#     for sampleA in namesList:
#         namesList = [a for a in namesList if a != sampleA]
#         for sampleB in namesList:
#             samplePairs.append((sampleA, sampleB))
#             outString += "\t" +f"{sampleA}-{sampleB}"
# print(outString)
#
# for chrom in outputDict:
#     for index in outputDict[chrom]:
#         if args.output_prefix is None:
#             currentSeqRecords = [SeqRecord(Seq(outputDict[chrom][index][a]), id=a) for a in outputDict[chrom][index]]
#             currentAlignment = MultipleSeqAlignment(records=currentSeqRecords)
#
#             if len(currentAlignment[0]) > 0:
#
#                 # Construct NJ tree
#                 calculator = DistanceCalculator('identity')
#                 constructor = DistanceTreeConstructor(calculator, 'nj')
#                 tree = constructor.build_tree(currentAlignment)
#
#                 # Print the tree
#                 outString = f"{chrom}\t{index}\t{len(currentAlignment[0])}\t{str(tree.format(fmt='newick')).strip()}"
#                 if args.distance:
#                     for samplePair in samplePairs:
#                         outString += f"\t" + f"{tree.distance(samplePair[0], samplePair[1])}"
#                 print(outString)
#                 # Calculate distances if needed
#
#         #else:
#         #    with open(f"{args.output_prefix}_{chrom}_{index}.fasta", "w") as outputFasta:
#         #        for sample in outputDict[chrom][index]:
#         #            outputFasta.write(f">{sample}\n{outputDict[chrom][index][sample]}\n")
