#!/bin/python
import os
import argparse
import csv
import numpy as np

def outputFormat(blocks):
    return [blocks[0][0], blocks[0][1], blocks[::-1][0][2]] +\
           [str(c) for c in np.sum([[int(b) for b in a[3:]] for a in blocks], axis=0)]


parser = argparse.ArgumentParser(description='Consolidate spectra profile to a larger window size. Output size must be'
                                             ' a multiple of Input size')
parser.add_argument('-i', dest='spectra_tsv', type=str, help='Input spectra tsv', required=True)
parser.add_argument('-m', dest='input_size', type=int, help='Input spectra window size', default=1000)
parser.add_argument('-n', dest='output_size', type=int, help='Output spectra window size', default=20000)
parser.add_argument('-o', dest='report_output', type=str, help='Output file name', default='consolidated_spectra.tsv')
args = parser.parse_args()

# make sure both file paths exist
if not os.path.exists(args.spectra_tsv):
    print("Cannot find spectra profile: "+args.spectra.tsv)
    exit()

# Adjust input width and output width to accommodate proper proportions

with open(args.spectra_tsv, "r") as inFile, open(args.report_output, "w") as outFile:
    inReader = csv.reader(inFile, delimiter='\t')
    header = next(inReader)
    outFile.write("\t".join(header)+"\n")
    currentBlock = [next(inReader)]
    observedWindowSize = int(currentBlock[0][2]) - int(currentBlock[0][1]) + 1
    targetFactor = args.output_size // args.input_size\
        if observedWindowSize == args.input_size\
        else args.output_size // observedWindowSize
    i = 0
    for line in inReader:
        if currentBlock[0][0] != line[0]:
            # output
            outFile.write("\t".join(outputFormat(currentBlock))+"\n")
            currentBlock = [line]
            i = 0
        elif i == targetFactor-1:
            i = 0
            # output
            outFile.write("\t".join(outputFormat(currentBlock))+"\n")
            currentBlock = [line]
        else:
            i += 1
            currentBlock.append(line)
