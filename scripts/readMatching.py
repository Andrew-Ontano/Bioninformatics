#!/usr/bin/env python3

# Tools for taking inputs of:
#    raw read name list from 'rawReadTools head'
#    *.noseq.gfa output from HiFiASM
#    *.paf output from ragtag scaffolding
# and outputting a tsv with the detailed location of each raw read

import logging
import argparse
import os
import time
from Bio import SeqIO
import pandas as pd
import re
from readpaf import parse_paf

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='readMatching')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-r', '--rawreads', dest='input_raw', type=str, help='Raw read list file', required=True)
parser.add_argument('-g', '--gfa', dest='input_gfa', type=str, help='GFA file', required=True)
parser.add_argument('-p', '--paf', dest='input_paf', type=str, help='PAF file', required=True)
parser.add_argument('-o', '--output', dest='output_table', type=str, help='Output table file', required=True)

args = parser.parse_args()
startTime = time.time()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_raw) or not os.path.exists(args.input_gfa) or not os.path.exists(args.input_paf):
    logging.error(f"Couldn't find one or more input files")
    exit()

read_list = []
with open(args.input_raw, 'r') as input_raw:
    read_list = [a.strip() for a in input_raw.readlines()]

paf_dataframe = pd.DataFrame()
with open(args.input_paf, 'r') as input_paf:
    paf_dataframe = parse_paf(input_paf, dataframe=True)

# read in the input gfa and parse out only (A)lignments
gfa_dataframe = pd.read_csv(args.input_gfa, sep='\t', names=['Type', 'ContigName', 'ContigStart', 'Strand', 'ReadName', 'ReadStart', 'ReadEnd', 'ReadID', 'Haplotype'])
gfa_dataframe = gfa_dataframe.loc[gfa_dataframe['Type'] == 'A']
with open(args.output_table, 'w') as output_table:
    output_table.write(f'Read\tContig\tContigStart\tContigEnd\tScaffold\tScaffoldStart\tScaffoldEnd\tStrand\n')
    for read_name in read_list:
        read_match = gfa_dataframe.loc[gfa_dataframe['ReadName'] == read_name]
        if len(read_match.index) == 0:
            # no matches
            output_table.write(f'{read_name}\t\t\t\t\t\t\t\n')
            pass
        else:
            for read_row in read_match.iterrows():
                gfa_match = paf_dataframe.loc[paf_dataframe['query_name'] == read_row[1]['ContigName']]
                if len(gfa_match.index) == 0:
                    output_table.write(f"{read_name}\t{read_row[1]['ContigName']}\t\t\t\t\t\t\n")
                else:
                    for gfa_row in gfa_match.iterrows():
                        output_table.write(f"{read_name}\t{read_row[1]['ContigName']}\t{read_row[1]['ContigStart']}\t{int(read_row[1]['ContigStart'])+int(read_row[1]['ReadEnd'])}\t{gfa_row[1]['target_name']}\t{gfa_row[1]['target_start']}\t{gfa_row[1]['target_end']}\t{'+' if read_row[1]['Strand'] == gfa_row[1]['strand'] else '-'}\n")
