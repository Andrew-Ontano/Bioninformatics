#!/bin/python

# given a fasta file, break each sequence into chunks of size S
import sys
import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='file', type=str, help='Input file name',required=True)
parser.add_argument('-o', dest='outfile', type=str, help='Output file name', default=False)
parser.add_argument('-s', dest='size', type=int, help='Size of chunks', default=1000)
args = parser.parse_args()

if not os.path.exists(args.file):
    print("File parameter is incorrect.")
    exit()

genome = SeqIO.to_dict(SeqIO.parse(args.file,'fasta'))

if args.outfile:
    with open(args.outfile,"w") as f:
        for sequence in genome:
            index = 0
            length = len(genome[sequence])
            while index < length:
                tail_index = index + args.size
                if tail_index > length:
                    tail_index = length
                f.write(">" + sequence + "_" + str(index) + ":" + str(tail_index - 1)+"\n")
                f.write(str(genome[sequence][index:tail_index].seq)+"\n")
                index += args.size
else:
    for sequence in genome:
        index = 0
        length = len(genome[sequence])
        while index < length:
            tail_index = index + args.size
            if tail_index > length:
                tail_index = length
            print(">"+sequence+"_"+str(index)+":"+str(tail_index-1))
            print(genome[sequence][index:tail_index].seq)
            index += args.size