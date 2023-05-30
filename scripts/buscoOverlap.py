#!/usr/bin/env python3

# read in all busco tables
# assess which are always present
# assess which are always missing

#for each busco, build a set that has an entry for each library?

import logging
import argparse
import os

logging.basicConfig(level=logging.INFO)
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest='inputs', help='Input busco tsvs, separated by spaces', nargs='*', required=True)
parser.add_argument('-o', '--output', dest='size', type=int, help='Genome size. Default results in raw size report only.')
args = parser.parse_args()

def tabulate(valueList):
    return [len([a for a in valueList if a == 'Complete']), len([a for a in valueList if a == 'Duplicated']), len([a for a in valueList if a == 'Fragmented']), len([a for a in valueList if a == 'Missing'])]


if not all([os.path.isfile(f) for f in args.inputs]):
    logging.error('Some or all of the input files do not exist.')
    exit()

buscosIndividual = {}
buscosCollated = {}

for file in args.inputs:
    fileName = file.replace('table_', '').replace('.tsv', '')
    tempDict = {}
    with open(file, 'r') as fileProcessing:
        for line in fileProcessing:
            if line[0] != '#':
                lineSplit = line.split('\t')
                lineName = lineSplit[0].strip()
                lineStatus = lineSplit[1].strip()
                tempDict[lineName] = lineStatus

                if lineName in buscosCollated:
                    buscosCollated[lineName][fileName] = lineStatus
                else:
                    buscosCollated[lineName] = {}
                    buscosCollated[lineName][fileName] = lineStatus

    buscosIndividual[fileName] = tempDict


libraryNames = list(buscosIndividual.keys())
tails = ['20', '30', '40', '50']
heads = ['m64310e_220613_200348', 'm64310e_220615_052921', 'm64310e_220616_145648', 'm64310e_220618_015229']
categories = ['Complete', 'Duplicated', 'Fragmented', 'Missing']
###### Tabulates the states of each busco in the *20, *30, *40, *50, and for each cell
print(f"busco,{','.join([f'{b}_{a}' for b in tails for a in categories])},{','.join([f'{b}_{a}' for b in heads for a in categories])}")
for busco in buscosCollated:
    tailDict = {}
    for tail in tails:
        tailDict[tail] = tabulate(list({key: value for key, value in buscosCollated[busco].items() if key.endswith(tail)}.values()))
    headDict = {}
    for head in heads:
        headDict[head] = tabulate(list({key: value for key, value in buscosCollated[busco].items() if key.startswith(head)}.values()))
    print(f"{busco},{str([tailDict[a] for a in tailDict]).replace('[', '').replace(']', '').replace(' ', '')},{str([headDict[a] for a in headDict]).replace('[', '').replace(']', '').replace(' ', '')}")

###### Writes the table to stdout
# print(f"busco,{','.join(buscosIndividual.keys())}")
# for busco in buscosCollated:
#     print(f"{busco},{','.join(buscosCollated[busco].values())}")

