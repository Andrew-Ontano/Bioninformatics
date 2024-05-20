#!/usr/bin/env python3

import sys
import pandas as pd
import os
import logging
import numpy as np

logging.basicConfig(level=logging.INFO)

if len(sys.argv) < 2:
    logging.warning(f'No file supplied')
    exit()
if not os.path.exists(sys.argv[1]):
    logging.warning(f'File {sys.argv[1]} does not exist.')
    exit()

with open(sys.argv[1]) as inFile:
    data = pd.read_csv(inFile, sep='\t')

libraries = data['Library'].unique()

groupedData = [data[data['Library'] == a] for a in libraries]

outData = pd.DataFrame(columns=data.columns)

for index in range(len(groupedData[0])):
    if 0 not in [datum.iloc[index, 4:68].sum() for datum in groupedData]:
        for datum in groupedData:
            outData.loc[len(outData.index)] = datum.iloc[index, 0:68]
outData.to_csv(f'gapless_{sys.argv[1]}', sep='\t')
#for datum in data.iterrows():
#    logging.info(datum[1][4:].sum())
