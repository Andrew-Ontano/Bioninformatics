#!/bin/python
import csv
import numpy as np
import scipy.stats as sp
import argparse
import os

def getFrequencies(tsvFile):
    with open(tsvFile, 'r') as tallyFile:
        inReaderTally = csv.reader(tallyFile, delimiter='\t')
        headers = next(inReaderTally)

        frequencyGroup = [next(inReaderTally)]
        if headers[0] != "Species":
            species = os.path.basename(tsvFile)
            frequencyGroup[0].insert(0, species)
        # pass the windows through a "baseline" filter, removing all datum that do not exceed these averages
        frequencies = {}
        for line in inReaderTally:
            if species:
                line.insert(0, species)
            if frequencyGroup[0][0] != line[0] or frequencyGroup[0][1] != line[1]:
                # need to process
                colSums = np.sum([[int(b) for b in a[::-1][0:64][::-1]] for a in frequencyGroup], axis=0)
                totalSums = sum(colSums)
                try:
                    frequencies[frequencyGroup[0][0] + "_" + frequencyGroup[0][1]] = [a / totalSums for a in colSums]
                except ZeroDivisionError:
                    frequencies[frequencyGroup[0][0] + "_" + frequencyGroup[0][1]] = [0.0] * 64
                frequencyGroup = [line]
            else:
                frequencyGroup.append(line)
        try:
            colSums = np.sum([[int(b) for b in a[::-1][0:64][::-1]] for a in frequencyGroup], axis=0)
            totalSums = sum(colSums)
            frequencies[frequencyGroup[0][0] + "_" + frequencyGroup[0][1]] = [a / totalSums for a in colSums]
        except ZeroDivisionError:
            frequencies[frequencyGroup[0][0] + "_" + frequencyGroup[0][1]] = [0.0] * 64
    return frequencies


parser = argparse.ArgumentParser(description='Consolidate spectra profile to a larger window size. Output size must be'
                                             ' a multiple of Input size')
parser.add_argument('-i', dest='spectra_tsv', type=str, help='Input spectra tsv', required=True)
parser.add_argument('-o', dest='spectra_output', type=str, help='Output spectra tsv', default='output.tsv')
parser.add_argument('-r', dest='weightedRemoval', action='store_true', help='Remove windows with spectra frequencies similar to the whole sequence', default=False)
parser.add_argument('-n', dest='weightedSubtraction', action='store_true', help='Normalize spectra frequencies for each window by the frequencies for the whole sequence', default=False)
parser.add_argument('-d', dest='distanceGraph', type=str, help='Distance from regions of interest calculated from supplied table', default=False)

args = parser.parse_args()

# make sure both file paths exist
if not os.path.exists(args.spectra_tsv):
    print("Cannot find spectra profile: "+args.spectra_tsv)
    exit()
if args.distanceGraph and not os.path.exists(args.distanceGraph):
    print("Cannot find distance table: "+args.distanceGraph)
    exit()

tsvFile = 'C:/Users/andre/Documents/triplet_analysis/spectra_30k_concat.tsv'
centromereTable = 'C:/Users/andre/Documents/triplet_analysis/group_centromeres.csv'
updatedTsvFile = 'C:/Users/andre/Documents/triplet_analysis/spectra_32k_concat_updated.tsv'
weightedRemoval = False
weightedSubtraction = False
distanceGraph = False

centroTable = []
if args.distanceGraph:
    with open(args.distanceGraph, "r") as inCentro:
        centroReader = csv.reader(inCentro, delimiter=",")
        next(centroReader)
        for row in centroReader:
            centroTable.append(row)

if args.weightedRemoval or args.weightedSubtraction:
    frequencies = getFrequencies(args.spectra_tsv)
else:
    frequencies = {}

with open(args.spectra_tsv, 'r') as inFile, open(args.spectra_output, 'w') as outFile:
    inReader = csv.reader(inFile, delimiter='\t')
    header = next(inReader)
    if header[0] != "Species":
        species = os.path.basename(args.spectra_tsv)
        header.insert(0, "Species")
    outFile.write("\t".join(header) + ("\tDistance\n" if args.distanceGraph else "\n"))
    for line in inReader:
        values = [int(a) for a in line[::-1][0:64][::-1]]
        rowSum = sum(values)
        newValues = [a/rowSum if rowSum != 0 else 0.0 for a in values]
        if args.weightedSubtraction:
            if species:
                line.insert(0, species)
            oldValues = [a for a in newValues]
            #newValues = [abs(oldValues[a]-frequencies[line[0]+"_"+line[1]][a]) for a in range(0, len(oldValues))]
            newValues = [abs(oldValues[a]-frequencies[line[0]+"_"+line[1]][a]) if oldValues[a] >= frequencies[line[0]+"_"+line[1]][a] else 0 for a in range(0, len(oldValues))]
        #newValues = []
        #for a in values:
        #    try:
        #        newValues.append(str(a/rowSum))
        #    except ZeroDivisionError:
        #        newValues.append("0")
        if args.distanceGraph:
            lineDistances = [(int(a[1]), int(a[2])) for a in centroTable if a[0] == line[1] and a[1] != ""]
            if len(lineDistances) > 0:
                lineDistance = "0"
                if lineDistances[0][0] > int(line[3]):
                    lineDistance = str(lineDistances[0][0] - int(line[3]))
                elif lineDistances[0][1] < int(line[2]):
                    lineDistance = str(int(line[2]) - lineDistances[0][1])
            else:
                lineDistance = "NA"

        outRow = "\t".join(line[0:4])
        if args.weightedRemoval:
            chiValues = [a for a in oldValues] if args.weightedSubtraction else [a for a in newValues]
            if sum(chiValues) != 0:
                chiObserved = [a for a in frequencies[line[0] + "_" + line[1]]]
                chiResult = sp.chisquare(chiValues, f_exp=chiObserved)[1]
                if chiResult <= 0.99:
                    continue

        outFile.write("\t".join(line[0:4]) + "\t" + "\t".join([str(a) for a in newValues]) + ("\t" + lineDistance + "\n" if args.distanceGraph else "\n"))

