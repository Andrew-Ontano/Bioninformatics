#!/usr/bin/env python3
import os
import argparse

parser = argparse.ArgumentParser(description='Apply labels to description of gff')
parser.add_argument('-i', dest='input_gff', type=str, help='Input gff file', required=True)
parser.add_argument('-d', dest='descriptions', type=str, help='Input descriptions file', required=True)
parser.add_argument('-o', dest='output_gff', type=str, help='Output gff file', default='new.gff')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
args = parser.parse_args()

if not os.path.exists(args.input_gff) or not os.path.exists(args.descriptions):
    print("Error: input file(s) not found")
    exit()

with open(args.output_gff, 'w') as output, open(args.input_gff, 'r') as gff, open(args.descriptions, 'r') as descriptions:
    descriptionDict = {}
    currentDescription = ""
    currentInfo = []
    for description in descriptions.readlines():
        if not description.startswith("#"):
            descriptionList = description.strip().split('\t')
            if currentDescription == descriptionList[0]:
                if descriptionList[8] not in currentInfo:
                    currentInfo.append(descriptionList[8])
            elif currentDescription != "":
                descriptionDict[currentDescription] = currentInfo
                currentInfo = []
                currentDescription = descriptionList[0]
            else:
                currentDescription = descriptionList[0]
        elif description.strip() == "##FASTA":
            descriptionDict[currentDescription] = currentInfo
            break
    for gffLine in gff.readlines():
        if gffLine.startswith("#"):
            output.write(gffLine)
        else:
            gffLineList = gffLine.strip().split("\t")
            currentAnnotation = gffLineList[8].split(";")
            if  gffLineList[8].startswith(currentAnnotation[0]) and currentAnnotation[0][3:] in descriptionDict:
                for annotation in descriptionDict[currentAnnotation[0][3:]]:
                    gffLineList[8] = annotation + "\n"
                    output.write("\t".join(gffLineList))
            else:
                output.write(gffLine)
