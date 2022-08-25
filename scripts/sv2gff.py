import argparse
import os
import logging
import csv

logging.basicConfig(level=logging.INFO)

# Tool for conversion of SVcaller outputs to gff. Designed for use in conjuction with Spectra, but is ammenable to any usage of GFF3
# Given the input SVcaller tsv: outputs a gff3 with the following column conversions:
# chrom ->  seqid
# ""    ->  source
# type  ->  type
# start ->  start
# end   ->  end
# refLength -> score
# "" -> strand
# "" -> phase
# "" -> attributes

parser = argparse.ArgumentParser(description='svcaller2gff file conversion')
parser.add_argument('-i', '--input', dest='input_sv', type=str, help='Input sv file', required=True)
parser.add_argument('-o', '--output', dest='output_gff', type=str, help='Output gff3 file', default='gff_output.gff')
parser.add_argument('-f', '--fasta-output', dest='fasta', action='store_true', help='Output fasta file of sv motifs', default=False)
parser.add_argument('-r', '--reference', dest='reference', action='store_true', help='Output reference positions instead of query positions', default=False)
args = parser.parse_args()

if not os.path.exists(args.input_sv):
    logging.error(f"Couldn't find input file '{args.input_sequence}'")
    exit()

#positions is a tuple of values to be referenced for line indices
if args.reference:
    positions = (5, 6, 7)
else:
    positions = (2, 3, 4)

with open(args.input_sv, 'r') as svFile:
    svContent = csv.reader(svFile, delimiter='\t')
    next(svContent)
    with open(args.output_gff, 'w') as gffFile:
        gffFile.write(f'##gff-version 3\n# file: {args.output_gff} derived from SVCaller file: {args.input_sv}\n')

        if args.fasta:
            with open(f'{args.output_gff}.fa', 'w') as fastaFile:
                for line in svContent:
                    gffFile.write(f'{line[0]}\tSVcaller\t{line[8]}\t{line[positions[0]]}\t{line[positions[1]]}\t{line[positions[2]]}\t+\t*\t*\n')
                    #gffFile.write([line[0], 'svcaller', line[8], line[2], line[3], line[4], '+', '*', '*'])
                    fastaFile.write(f'>{line[0]}_{line[positions[positions[1]]]}-{line[positions[2]]}_{line[8]}\n')
                    fastaFile.write(f'{line[10]}\n')

        else:
            for line in svContent:
                gffFile.write(f'{line[0]}\tSVcaller\t{line[8]}\t{line[positions[0]]}\t{line[positions[1]]}\t{line[positions[2]]}\t+\t*\t*\n')
