#!/bin/python3
import sys
from math import log10, ceil

if len(sys.argv)<2:
    print("No file argument given")
    sys.exit()
else:
    inFile = sys.argv[1]
    outFile = inFile+".fasta"

with open(inFile) as f:
    lines = f.readlines()

uid_size = ceil(log10(len(lines)))+1
uid_start = 1

fileWrite = open(outFile,"w")
for line in lines:
    fileWrite.write(">repeat_"+str(uid_start).zfill(uid_size)+"\n"+line)
    uid_start += 1
fileWrite.close()
