#!/bin/python

# important future plans
# two inputs: gfa and rawreads file list;
# populate seqBases with raw read sequences instead of corresponding contig
import sys

if len(sys.argv) < 2:
    print("GFA file path needed")
    quit()

fileIn = sys.argv[1]
fileOut = sys.argv[1]+".sam"
with open(fileIn, 'r') as f:
    contents = f.readlines()

seqDict = {}
with open(fileOut, 'w') as f:
    f.write("@HD VN:1.6 SO:coordinate\n")
    for seq in [a for a in contents if a[0] == "S"]:
        seqList = seq.split('\t')
        seqDict[seqList[1]] = seqList[2]
        lineOut = "@SQ SN:" + seqList[1] + " LN:" + seqList[3][5:] + "\n"
        f.write(lineOut)
    for seq in [a for a in contents if a[0] == "A"]:
        seqList = seq.split('\t')
        flag = "0" if seqList[3] == "+" else "99"
        seqBases = seqDict[seqList[1]][int(seqList[2]):int(seqList[2])+int(seqList[6])]
        lineOut = "\t".join([seqList[4], flag, seqList[1], seqList[2], "255", "*", "*", "0", str(int(seqList[6]) + 1), seqBases, "*"])+"\n"
        f.write(lineOut)
