#!/usr/bin/env python3
import argparse
from Bio import SeqIO,Seq

# Given an input fasta file of repeat elements, reports back a new fasta file with only the unique permutations of all element
def rc(seq):
    seq = Seq.Seq(seq)
    return str(seq.reverse_complement())


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parser.add_argument('-o', '--output', dest='output_sequence', type=str, help='Output sequence file', required=True)
args = parser.parse_args()

newSeqs = []
sequenceDict = SeqIO.to_dict(SeqIO.parse(args.input_sequence, format='fasta'))
for sequence in sequenceDict:
    sequenceLength = len(sequenceDict[sequence].seq)
    rotations = []
    for a in range(sequenceLength):
        testSeq = f'{sequenceDict[sequence].seq[a:]}{sequenceDict[sequence].seq[0:a]}'
        rotations.append(testSeq)
        rotations.append(rc(testSeq))
    if any([a in newSeqs for a in rotations]):
        print()
    else:
        newSeqs.append(str(sequenceDict[sequence].seq))
    print(rotations)
print(newSeqs)

