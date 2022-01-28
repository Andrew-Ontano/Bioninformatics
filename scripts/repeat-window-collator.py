#!/bin/python
import sys

if len(sys.argv)<2:
    print("tsv file path needed")
    quit()

fileIn = sys.argv[1]
fileOut = sys.argv[1]+".sam"
with open(fileIn,'r') as f:
    contents = f.readlines()

# define the row headers for safekeeping, then remove the header row
headerRow = contents[0].strip().split('\t')
headers = ["Scaffold","Length"]+headerRow[3:len(headerRow)]
contents.pop(0)

entries = {}

for content in contents:
    row =  content.strip().split('\t')
    if row[0] in entries:
        # add to row
        entries[row[0]][0] = int(row[2])
        entries[row[0]][1:] = [a + int(b) for a, b in zip(entries[row[0]][1:], row[3:])]
    else:
        # set the scaffold in the dictionary to the end seq plus each repeat count
        entries[row[0]] = [int(a) for a in row[2:]]

print('\t'.join(headers))
for entry in entries:
    print(entry+'\t'+'\t'.join(str(a) for a in entries[entry]))