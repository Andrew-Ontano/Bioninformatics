#!/usr/bin/env python3

# Given a txt file of read lengths (one per line) and a reference genome size, calculate recommend coverage cutoffs

import logging
import argparse
import os

# setup for cl args and logging
logging.basicConfig(level=logging.INFO)
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='inFile', type=str, help='Input file name.', required=True)
parser.add_argument('-s', dest='size', type=int, help='Genome size. Default results in raw size report only.')
args = parser.parse_args()

# coverage calculation
def getCoverage(lengthList, terminus, referenceGenomeSize):
    return sum(lengths[0:terminus])/referenceGenomeSize

# check for file
if not os.path.exists(args.inFile):
    logging.warning("File parameter is incorrect.")
    exit()

# start list and populate from file
lengths = list()
with open(args.inFile, 'r') as inFile:
    for line in inFile:
        lengths.append(int(line))

totalLength = sum(lengths)
if isinstance(args.size, int):
    # sort lengths to downsample by length
    lengths.sort(reverse=True)
    #get size estimates as seed points
    seeders = [.001, .003, .005, .007, .009, .01, .02, .03, .04, .05, .06, .07, .08, .09, .1,
               .111, .122, .133, .144, .155, .166, .177, .188, .2,
               .22, .24, .26, .28, .3, .32, .34, .36, .38, .4,
               .43, .46, .49, .52, .55, .58, .61, .64, .67, .7,
               .74, .78, .82, .86, .9, .94, .98, .99, .995, .999, .9995, 1]
    estimators = [int(round(a * len(lengths))) for a in seeders]
    placements = [getCoverage(lengths, a, args.size) for a in estimators]

    # get max coverage for cutting off impossible coverages
    maxCoverage = totalLength / args.size
    targets = [a for a in [20, 30, 40, 50] if a <= maxCoverage]

    # populate seed points
    seedPoints = {}
    for target in targets:
        seedPoints[target] = [estimators[placements.index(a)] for a in placements if a >= target][0]

    logging.info(f'Starting targets are {[lengths[seedPoints[a]] for a in seedPoints]} for targets {targets}')

    for seed in seedPoints:
        previousIndex = seedPoints[seed]
        while True:
            coverage = getCoverage(lengths, previousIndex-1, args.size)
            if coverage < seed:
                logging.info(f'Length {lengths[previousIndex]} is the size for {seed}.')
                logging.info(f'Actual coverage: {getCoverage(lengths, previousIndex, args.size)}')
                break
            else:
                previousIndex -= 1
else:
    # just output size
    logging.info(f'Raw genomic size: {totalLength}')
