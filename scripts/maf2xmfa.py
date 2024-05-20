#!/usr/bin/env python3

from Bio import AlignIO
import argparse
import os
import logging
import textwrap
import io

def prepareSequenceInfo(seqRecord):
    if seqRecord.annotations['strand'] == 1:
        return (
            seqRecord.annotations['start'] + 1,
            seqRecord.annotations['start'] + seqRecord.annotations['size'],
            '+'
        )
    else:
        return (
            seqRecord.annotations['srcSize'] - seqRecord.annotations['start'] - seqRecord.annotations['size'] + 1,
            seqRecord.annotations['srcSize'] - seqRecord.annotations['start'],
            '-'
        )


# Set up logger
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

parser = argparse.ArgumentParser(description='maf2xmfa: a tool for converting maf to xmfa.')
parser.add_argument('-i', '--input', dest='input_maf', type=str, help='Input MAF file', required=True)
parser.add_argument('-o', '--output', dest='output_xmfa', type=str, help='Output XMFA file', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-g', '--gappy', dest='gappy', action='store_true', help='Include gappy blocks', default=False)
args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_maf):
    logger.error(f"Couldn't find input MAF '{args.input_maf}'")
    exit()

collection = {}
blackList = {'Anc0', '_MINIGRAPH_'}
alignments = AlignIO.parse(args.input_maf, 'maf')

with open(args.output_xmfa, 'wb') as output_file, io.BufferedWriter(output_file) as buffered_output:
    logger.info(f"Processing file '{args.input_maf}'")
    # Collect alignment information and create name list
    for alignment in alignments:
        for align in alignment:
            name = align.id
            if name not in collection.keys() and name.split('.')[0] not in blackList:
                collection[name] = len(collection) + 1
    logger.info(f"Alignment names collected from '{args.input_maf}'")

    # Write dummy block and Mauve file headers
    output_string = "#FormatVersion Mauve1\n"
    second_string = ""
    for index, collectionName in enumerate(collection):
        output_string += f"#Sequence{index+1}Entry {collectionName}\n"
        second_string += f"> {collection[collectionName]}:0-0 + {collectionName}\n-\n"
    second_string += "=\n"
    output_string += second_string
    buffered_output.write(output_string.encode('utf-8'))

    logger.info(f"XMFA headers written to '{args.output_xmfa}'")

    alignments = AlignIO.parse(args.input_maf, 'maf')
    for alignment in alignments:
        names = {a.id for a in alignment if a.id.split('.')[0] not in blackList}

        # Write test validates whether any sequences were in an LCB, or all were blacklisted/missing
        written = False

        # If a sequence has already been reported in a block, the looped sequences should be added to a new block alone
        extraAlignments = []
        output_string = ""

        for collectionName in collection:
            if collectionName in names:
                written = True
                row = [a for a in alignment if a.id == collectionName]
                values = prepareSequenceInfo(row[0])
                output_string += f"> {collection[collectionName]}:{values[0]}-{values[1]} {values[2]} {collectionName}\n{textwrap.fill(str(row[0].seq), break_long_words=True, break_on_hyphens=False, width=80)}\n"

                if len(row) > 1:
                    extraAlignments.extend(row[1:])
            elif args.gappy:
                written = True
                # write a blank row
                output_string += f"> {collection[collectionName]}:0-0 + {collectionName}\n{textwrap.fill('-' * len(alignment[0]), break_long_words=True, break_on_hyphens=False, width=80)}\n"

        # Only write the footer if an LCB needs to be written
        if written:
            output_string += '=\n'

        if len(extraAlignments) > 0:
            for extraEntry in extraAlignments:
                values = prepareSequenceInfo(extraEntry)
                output_string += f"> {collection[extraEntry.id]}:{values[0]}-{values[1]} {values[2]} {extraEntry.id}\n{textwrap.fill(str(extraEntry.seq), break_long_words=True, break_on_hyphens=False, width=80)}\n=\n"

        if len(output_string) > 0:
            buffered_output.write(output_string.encode('utf-8'))

    logger.info(f"File written to '{args.output_xmfa}', exiting.")
