#!/usr/bin/env python3
import argparse
import os
from Bio import SeqIO, SeqRecord, pairwise2
from Bio.pairwise2 import format_alignment
import logging
from collections import Counter
from itertools import combinations
import numpy as np

def alignment_dict_to_matrix(alignment_dict):
    seq_names = {name for pair in alignment_dict.keys() for name in pair}
    score_matrix = [[0] * len(seq_names) for i in range(len(seq_names))]
    for i, name1 in enumerate(seq_names):
        for j, name2 in enumerate(seq_names):
            if name1 != name2:
                # score_matrix[i][j] = alignment_dict.get((name1, name2), alignment_dict.get((name2, name1), 0)).score
                score_matrix[i][j] = alignment_dict.get((name1, name2), alignment_dict.get((name2, name1), 0))
    score_matrix = np.array(score_matrix)
    return score_matrix

def write_score_matrix_to_tsv(score_matrix, seq_names, file_path):
    with open(file_path, 'w') as f:
        # Write header row with column labels
        gene_names = [gene.name for gene in genome]
        f.write("\t" + "\t".join(gene_names) + "\n")
        # Write data rows with row labels and scores
        for i in range(len(score_matrix)):
            row_string = gene_names[i] + "\t" + "\t".join([str(score) for score in score_matrix[i]]) + "\n"
            f.write(row_string)



def jaccard_similarity(sequenceA, sequenceB):
    intersection = len(sequenceA.intersection(sequenceB))
    union = len(sequenceA.union(sequenceB))
    return intersection / union

# setup for cl args and logging
logging.basicConfig(level=logging.INFO)
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='file', type=str, help='Input file name.', required=True)
parser.add_argument('-o', dest='out_file', type=str, help='Output file name.', default=None)
parser.add_argument('-m', dest='mer_size', type=int, help='Desired k-mer size', default=3)
args = parser.parse_args()

# check for file
if not os.path.exists(args.file):
    logging.error("File parameter is incorrect.")
    exit()

genome = []
for record in SeqIO.parse(args.file, "fasta"):
    genome.append(record)

sequence_similarity = {}

for pair in combinations(genome, 2):
    #test 5-5 pairing, then test 5-3 pairing
    alignment_forward = pairwise2.align.globalxx(pair[0].seq, pair[1].seq)[0]
    alignment_reverse = pairwise2.align.globalxx(pair[0].seq, pair[1].seq.reverse_complement())[0]
    if alignment_forward.score > alignment_reverse.score:
        sequence_similarity[(pair[0].name, pair[1].name)] = alignment_forward.score
    else:
        sequence_similarity[(pair[0].name, pair[1].name)] = alignment_reverse.score

if args.out_file:
    write_score_matrix_to_tsv(alignment_dict_to_matrix(sequence_similarity), genome, args.out_file)
else:
    print(alignment_dict_to_matrix(sequence_similarity))
