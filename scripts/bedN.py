#!/usr/bin/env python3

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='vcfFilter: a tool for processing vcf results by genotype information')
parser.add_argument('-i', '--input', dest='input_fasta', type=str, help='Input FASTA file', required=True)
parser.add_argument('-o', '--output', dest='output_bed', type=str, help='Output BED file', required=True)
args = parser.parse_args()

# Function to write a BED record
def write_bed_record(file, chrom, start, end, name=".", score="."):
    file.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\n")


# Open the input FASTA file and output BED file
with open(args.output_bed, "w") as bed_file:

    for record in SeqIO.parse(args.input_fasta, "fasta"):
        chrom = record.id  # Chromosome or contig name
        seq = str(record.seq)
        mask_start = None
        mask_end = None
        for i, base in enumerate(seq):
            if base == "N":
                # "N" indicates masked regions
                if mask_start is None:
                    mask_start = i
                mask_end = i
            elif mask_start is not None:

                write_bed_record(bed_file, chrom, mask_start, mask_end + 1)
                mask_start = None
                mask_end = None

        # Check if a masked region ends at the end of the sequence
        if mask_start is not None:
            write_bed_record(bed_file, chrom, mask_start, mask_end + 1)
