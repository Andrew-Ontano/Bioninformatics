#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description='Spectra genetic profiling. Counting, processing, and visualization of 3-mers')
parser.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
parser.add_argument('-o', '--output', '--output', dest='output', type=str, help='Output spectra tsv', default=False)
parser.add_argument('-r', '--weighted-filter', dest='weighted_filter', action='store_true',
                    help='Produce two additional outputs that have outlier windows and normal windows', default=False)
parser.add_argument('-n', '--weighted-norm', dest='weighted_normalization', action='store_true',
                    help='Normalize spectra frequencies for each window by the frequencies for the whole sequence',
                    default=False)
parser.add_argument('-f', '--freq', dest='frequencies', action='store_true',
                    help='Mark this is Spectra data is already in frequencies', default=False)
parser.add_argument('-c', '--convert', dest='convert', action='store_true',
                    help='Convert between counts and frequencies', default=False)
parser.add_argument('-s', '--window-resize', dest='resize_window', type=int,
                    help='Resize windows to summarize N for every 1 window')
parser.add_argument('-p', '--print', dest='print', action='store_true', help='Print global frequencies', default=False)
parser.add_argument('-y', '--simplify', dest='simplify', action='store_true',
                    help='Simplify forward and reverse-complement counts per window', default=False)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='verbose mode', default=False)

args = parser.parse_args()

spectra = pd.read_csv(args.input_tsv, delimiter='\t')
targetFactor = int(args.resize_window)
if targetFactor > 1:
    spectraGroups = spectra.groupby(['Sequence'])
    spectra = pd.DataFrame()
    for group in spectraGroups:
        for index in range(0, len(group[1]), targetFactor):
            subset = group[1].iloc[index:index + targetFactor]

            newSubset = pd.DataFrame({'Sequence': subset.iloc[0,0], 'Start': subset['Position'].min(), 'End': subset['Position'].max(), 'Cov': subset['Cov'].mean()}, index=[0])

            spectra = pd.concat([spectra, newSubset], ignore_index=True)
spectra.reindex()
spectra.to_csv(args.output, sep='\t', index=False)


