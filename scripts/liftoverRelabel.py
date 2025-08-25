#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from Bio import SeqIO, Seq

parser = argparse.ArgumentParser(description='Apply labels to description of gff')
parser.add_argument('-i', dest='input_vcf', type=str, help='Input vcf file', required=True)
parser.add_argument('-f', dest='reference_fasta', type=str, help='Input reference fasta file', required=True)
parser.add_argument('-c', dest='coordinates', type=str, help='Input coordinate file', required=True)
parser.add_argument('-o', dest='output_vcf', type=str, help='Output vcf file', default='new.vcf')
parser.add_argument('-e', dest='output_error', type=str, help='Output error file', default='new.error')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
args = parser.parse_args()

if not os.path.exists(args.input_vcf) or not os.path.exists(args.coordinates) or not os.path.exists(args.reference_fasta):
    print("Error: input file(s) not found")
    exit()

referenceGenome = SeqIO.to_dict(SeqIO.parse(args.reference_fasta, format="fasta"))
with open(args.output_vcf, 'w') as output, open(args.output_error, 'w') as error, open(args.input_vcf, 'r') as vcf, open(args.coordinates, 'r') as coordinates:
    coords = pd.read_csv(coordinates, dtype=str, sep='\t', header=None, comment="#")


    for line in vcf.readlines():
        if line.startswith("#"):
            output.write(line)
        else:
            temp = line.strip().split('\t')
            tempDF = coords[coords[11] == temp[1]]

            if len(tempDF) == 0:
                error.write(line)
            else:
                # get lookup position from reference
                referenceBase = referenceGenome[tempDF.iloc[0,13]].seq[int(tempDF.iloc[0,15])+1]
                # RC variants if necessary
                temp[0] = tempDF.iloc[0, 13]
                temp[1] = tempDF.iloc[0, 15]

                if tempDF.iloc[0,8] == "+-" and str(Seq.Seq(referenceBase).reverse_complement()).lower() == temp[3].lower():
                    temp[3] = referenceBase
                    temp[4] = str(Seq.Seq(temp[4]).reverse_complement())
                elif tempDF.iloc[0,8] == "++" and referenceBase.lower() != temp[3].lower():
                    if referenceBase == temp[4]:
                        replaceTemp = []
                        for a in temp[8:len(temp)]:
                            tempField = a.split(":")
                            tempField[0].replace("1","q").replace("1","0").replace("q","1")
                            replaceTemp.append(":".join(tempField))
                        #replaceTemp = [a.split(":")[0].replace("0","q").replace("1","0").replace("q","1") for a in temp[8:len(temp)]]
                        temp[8:len(temp)] = replaceTemp
                        temp[4] = temp[3]
                        temp[3] = referenceBase
                    else:
                        replaceTemp = []
                        for a in temp[8:len(temp)]:
                            tempField = a.split(":")
                            tempField[0].replace("1", "2").replace("0", "1")
                            replaceTemp.append(":".join(tempField))
                        #replaceTemp = [a.split(":")[0].replace("0", "2") for a in temp[8:len(temp)]]
                        temp[8:len(temp)] = replaceTemp
                        temp[4] = f"{temp[3]},{temp[4]}"
                        temp[3] = referenceBase
                elif tempDF.iloc[0,8] == "+-" and str(Seq.Seq(referenceBase).reverse_complement()).lower() != temp[3].lower():
                    if str(Seq.Seq(referenceBase).reverse_complement()) == temp[4]:
                        replaceTemp = []
                        for a in temp[8:len(temp)]:
                            tempField = a.split(":")
                            tempField[0].replace("0","q").replace("1","0").replace("q","1")
                            replaceTemp.append(":".join(tempField))
                        #replaceTemp = [a.split(":")[0].replace("0","q").replace("1","0").replace("q","1") for a in temp[8:len(temp)]]
                        temp[8:len(temp)] = replaceTemp
                        temp[4] = str(Seq.Seq(temp[3]).reverse_complement())
                        temp[3] = referenceBase
                    else:
                        replaceTemp = []
                        for a in temp[8:len(temp)]:
                            tempField = a.split(":")
                            tempField[0].replace("1", "2").replace("0", "1")
                            replaceTemp.append(":".join(tempField))
                        #replaceTemp = [a.split(":")[0].replace("0", "2") for a in temp[8:len(temp)]]
                        temp[8:len(temp)] = replaceTemp
                        temp[4] = f"{str(Seq.Seq(temp[3]).reverse_complement())},{str(Seq.Seq(temp[4]).reverse_complement())}"
                        temp[3] = referenceBase
                output.write("\t".join(temp)+"\n")




