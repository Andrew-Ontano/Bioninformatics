#!/usr/bin/env python3

# tableProcess is a generalized script for wrapping a human-readable table with pandas commands. This will allow you to filter the table by row operations

import logging
import argparse
import os
import time
from Bio import SeqIO, Seq
import re
import pandas as pd

# Set up logger
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='tableProcess: a tool for processing tables using pandas')

parser.add_argument('-i', '--input', dest='input_table', type=str, help='Input table file', required=True)
parser.add_argument('-i2', '--input2', dest='input_table2', type=str, help='Input table file', required=False)
parser.add_argument('-o', '--output', dest='output_table', type=str, help='Output table file', required=True)
parser.add_argument('-c', '--command', dest='pandas_command', type=str, help='Command to give pandas', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-s', '--separator', dest='table_separator', type=str, help='Table separator value', default='\t')
args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_table):
    logging.error(f"Couldn't find input file '{args.input_table}'")
    exit()

table = pd.read_csv(args.input_table, delimiter=args.table_separator)

# define table with vcf names
#table = pd.read_csv(args.input_table, delimiter=args.table_separator, names=['Chrom', 'Pos', 'ID', 'Ref', 'Alt', 'Qual', 'Filter', 'Info', 'Format', 'VmandF', 'Vpenn', 'Vvelu'])

# filter trf tsv based on region length
#outputTable = table.loc[table['region_sequence'].str.len() > 5000]

# filter vcf for variants longer than 30 in either ref or alt
#outputTable = table.loc[table['Ref'].str.len() > 29]
#outputTable = table.loc[table['Alt'].str.len() > 29]

#outputTable.to_csv(args.output_table, sep=args.table_separator, index=False)

# merge two tables
##df1 = pd.read_csv(args.input_table, delimiter=args.table_separator, header=None)
##df2 = pd.read_csv(args.input_table2, delimiter=args.table_separator, header=None)
##df3 = pd.concat([df1, df2], ignore_index=True).drop_duplicates()
##df3.to_csv(args.output_table, sep=args.table_separator, index=False)

# filter results that are consistent with the outgroup
##outputTable = table.loc[table['VmandF'] != '0']
##outputTable = outputTable.loc[outputTable['VmandF'] != '.']

# Add an amount specificed with 'command' to each start/stop value
table['start'] = table['start'] + int(args.pandas_command)
table['end'] = table['end'] + int(args.pandas_command)

table.to_csv(args.output_table, sep=args.table_separator, index=False)

# general output
#outputTable.to_csv(args.output_table, sep=args.table_separator, index=False)
