#!/bin/python3
import os
import argparse
from Bio import SeqIO
import csv
from math import ceil, floor

def getFrequencies(values):
    if type(values[0]) == list:
        blocks = [sum([float(j) for j in i]) for i in zip(*values)]
        return getFrequencies(blocks)
    else:
        numerics = [float(a) for a in values]
        total = sum(numerics)
        try:
            return [b/total for b in numerics]
        except ZeroDivisionError:
            return [0] * len(numerics)

def processBlocks(values):
    if len(values) == 1:
        return values
    entries = [int(a) for a in list(values.keys())]
    entries.sort()
    current_block_start = str(entries[0])
    initial_frequencies = getFrequencies(values[current_block_start])
    blocks = {}
    current_block = [values[str(current_block_start)]]

    for i in range(1, len(entries)):
        entry_key = str(entries[i])
        current_frequencies = getFrequencies(values[entry_key])
        frequency_delta = sum(
            [abs(initial_frequencies[a] - current_frequencies[a]) for a in range(len(current_frequencies))])
        if frequency_delta <= float(args.similarity):
            current_block.append(values[entry_key])
        else:
            blocks[current_block_start] = getFrequencies(current_block)
            current_block_start = entry_key
            current_block = [values[entry_key]]
        initial_frequencies = current_frequencies
    blocks[current_block_start] = getFrequencies(current_block)
    return blocks


def percentChange(a, b):
    if a < b:
        x, y = a, b
    else:
        x, y = b, a
    try:
        return abs(x - y) / y
    except ZeroDivisionError:
        return 0


parser = argparse.ArgumentParser(description='Analyze a spectra 3-mer output for congruent areas')
parser.add_argument('-i', dest='input_spectra', type=str, help='Input spectra tsv', required=True)
parser.add_argument('-o', dest='report_output', type=str, help='Output file name', default='spectra_block_report.tsv')
parser.add_argument('-p', dest='proportion', type=str, help='Proportional similarity  threshold for spectra frequencies', default='.1')
parser.add_argument('-m', dest='minimum_block_size', type=str, help='Minimum length of a block to be kept', default='12000')
parser.add_argument('-q', dest='quiet', action='store_false', help='Print info on the run', default=True)
parser.add_argument('-s', dest='similarity', type=str, help='Similarity threshold', default='0.1')
parser.add_argument('-r', dest='refine', action='store_true', help='block refinement passes', default=False)

args = parser.parse_args()

if not os.path.exists(args.input_spectra):
    print("Error: input file not found")
    exit()


"""
    Read through tsv, creating a dictionary of each sequence containing a dictionary of spectra frequencies keyed by the
    start position of the window
"""

spectra_data = {}
window_size = 0

seq_lengths = {}

with open(args.input_spectra) as f:
    spectra_tsv = csv.reader(f, delimiter='\t')
    header = next(spectra_tsv)
    for row in spectra_tsv:
        if not window_size:
            window_size = int(row[2]) - int(row[1]) + 1
        row_label = row[0]
        if row_label in spectra_data:
            spectra_data[row_label][row[1]] = row[3:len(header)]
            seq_lengths[row_label] = row[2]
        else:
            spectra_data[row_label] = {}
            spectra_data[row_label][row[1]] = row[3:len(header)]

"""
    Process each sequence separately since observations are independent
"""

# Read through spectra file, and process each scaffold separately
# Find potential breakpoints
# Current code commented out
with open(args.report_output, 'w') as f:
    f.write("Library\tSequence\tPositions\t"+"\t".join(header[3:len(header)])+"\n")
    for key in spectra_data:
        if args.quiet:
            print("Processing " + key + "...")
        values = processBlocks(spectra_data[key])
        if args.refine:
            currentLength = len(values)
            values = processBlocks(values)
            while len(values) > 1:
                if args.quiet:
                    print("\tRefining...")
                values = processBlocks(values)
                if len(values) >= currentLength:
                    break
                else:
                    currentLength = len(values)

        keys = list(values.keys())
        for i in range(0, len(keys)-1):
            block_positions = str(keys[i]) + "-" + str(int(keys[i+1])-1)
            if args.quiet:
                print("\tPositions: " + block_positions)
            f.write(args.input_spectra + "\t" + key + "\t" + block_positions + "\t" + "\t".join(
                [str(a) for a in values[keys[i]]]) + "\n")
        block_positions = str(keys[::-1][0]) + "-" + seq_lengths[key]
        if args.quiet:
            print("\tPositions: " + block_positions)
        f.write(args.input_spectra + "\t" + key + "\t" + block_positions + "\t" + "\t".join(
            [str(a) for a in values[keys[::-1][0]]]) + "\n"
        )

"""
with open(args.report_output, 'w') as f:
    f.write("Library\tSequence\tPositions\t"+"\t".join(header[3:len(header)])+"\n")

    for key in spectra_data:
        entries = [int(a) for a in list(spectra_data[key].keys())]
        entries.sort()

        entry = entries[0]

        if args.quiet:
            print("Processing " + key + "...")
        initial_frequencies = getFrequencies(spectra_data[key][str(entry)])
        current_block_length = 1
        current_block_start = entry

        for entry in entries[1:]:
            entry_frequencies = getFrequencies(spectra_data[key][str(entry)])
            frequency_deltas = [percentChange(initial_frequencies[a], entry_frequencies[a]) for a in range(0, len(initial_frequencies))]
            if sum([a <= float(args.proportion) for a in frequency_deltas]) >= len(entry_frequencies)*float(args.similarity):
                initial_frequencies = [((initial_frequencies[a] * current_block_length) + entry_frequencies[a]) / (current_block_length+1) for a in range(0, len(initial_frequencies))]
                current_block_length += 1
            else:
                if (current_block_length * window_size) >= int(args.minimum_block_size):
                    block_positions = str(current_block_start) + "-" + str(current_block_start - 1 + current_block_length * window_size)
                    if args.quiet:
                        print("Positions: " + block_positions)
                    f.write(args.input_spectra+"\t"+key+"\t"+block_positions+"\t"+"\t".join([str(a) for a in initial_frequencies])+"\n")
                initial_frequencies = entry_frequencies
                current_block_length = 1
                current_block_start = entry
        if (current_block_length * window_size) >= int(args.minimum_block_size):
            block_positions = str(current_block_start) + "-" + seq_lengths[key]
            if args.quiet:
                print("Positions: " + block_positions)
            f.write(args.input_spectra + "\t" + key + "\t" + block_positions + "\t" + "\t".join([str(a) for a in initial_frequencies])+"\n")
"""
