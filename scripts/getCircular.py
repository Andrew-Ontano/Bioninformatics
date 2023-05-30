from Bio import Seq
import sys

def rc(sequence):
    sequence = Seq.Seq(sequence)
    return str(sequence.reverse_complement())

def rotate(string, i=1):
    return f'{string[i:]}{string[0:i]}'


with open(sys.argv[1], 'r') as seqListFile:
    seqList = seqListFile.readlines()
cleanList = []
for seq in seqList:
    seq = seq.strip()
    rotations = [rotate(seq, a) for a in range(len(seq))] + [rotate(rc(seq), a) for a in range(len(seq))]
    breaker = True
    for check in rotations:
        if check in cleanList:
            breaker = False
            break
    if breaker:
        cleanList.append(seq)

print(','.join(cleanList))
