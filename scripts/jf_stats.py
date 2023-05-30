# Given an input file of kmer counts from jellyfish2, outputs median kmer count

import sys
from statistics import median
from Bio import SeqIO
df = sequenceDict = SeqIO.to_dict(SeqIO.parse(sys.argv[1]))

print(df)


#with open(sys.argv[1], 'r') as file:
#    fileContents = file.readlines()
#    merDict = {}
#    for a in range(0, len(fileContents), 2):
#        merDict[fileContents[a+1].strip()] = int(fileContents[a][1:].strip())

#print(sorted(merDict.values())[len(merDict)//2])
