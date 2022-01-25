import sys
from Bio import SeqIO
import argparse

def n50(sequenceLengths):
    workLengths = sorted(sequenceLengths,reverse=True)
    total = sum(sequenceLengths)
    tallyLength=0
    tallyNum=1
    for seq in workLengths:
        tallyLength=tallyLength+seq
        if tallyLength >= total/2:
            return(tallyNum)
        else:
            tallyNum=tallyNum+1
 
parser = argparse.ArgumentParser()
parser.add_argument('-in', dest='file', type=str, help='Input file name')
parser.add_argument('-type', dest='file_type', type=str, help='File type; default = f(asta); other options: currently unavailable', default="f")
parser.add_argument('-data', dest='data_type', type=str, help='Data type; default = s(equence); other options: a(lignment)', default="s")
parser.add_argument('-seq', dest='sequence_type', type=str, help='Sequence type; default = n(ucleotide); other options: p(eptide)', default="n")
parser.add_argument('-n50', dest='calculate_n50', action='store_true', help='Returns the N50 sequence length', default=False)
parser.add_argument('-bf', dest='calculate_base_freqs', action='store_true', help='Returns the frequency each of base', default=False)
parser.add_argument('-mpsi', dest='calculate_pairwise_identity', action='store_true', help='Returns the mean pairwise identity', default=False)
parser.add_argument('-quiet', dest='quiet_output', action='store_true', help='Omits basic sequence info', default=False)
args = parser.parse_args()
if not args.file:
    print("Please specify an input file with -in file")
    exit()

file_type = {"f":"fasta"}
sequence_type = {"n":"acgt","p":"arndceqghilkmfpstwyv"}

sequenceDict = SeqIO.to_dict(SeqIO.parse(args.file,file_type[args.file_type]))
if args.calculate_base_freqs or args.calculate_pairwise_identity:
    baseSums = dict.fromkeys(list(sequence_type[args.sequence_type]),0)
    sequenceDictKeys = list(sequenceDict.keys())
    mpsi = []
    for seq_record in sequenceDictKeys:
        if args.calculate_base_freqs:
            for base in baseSums:
                baseSums[base] += sum([a==base for a in sequenceDict[seq_record].lower()])
        if args.calculate_pairwise_identity and args.data_type=="a":
            for second_sequence in sequenceDictKeys[sequenceDictKeys.index(seq_record)+1:]:
                scoreSum = sum([1 for a in zip(enumerate(sequenceDict[seq_record].seq),enumerate(sequenceDict[second_sequence].seq)) if len(set(a))==1 and '-' not in str(a)])
                scoreLen = len(set([i for i, ltr in enumerate(sequenceDict[seq_record].seq) if ltr != "-"]).intersection(set([i for i, ltr in enumerate(sequenceDict[second_sequence].seq) if ltr != "-"])))
                if scoreLen != 0:
                    score = scoreSum/scoreLen
                    mpsi.append(score)


total = sum([len(sequenceDict[seq].seq) for seq in sequenceDict])

if args.quiet_output:
    outString = ""
else:
    outString = "Sequences: "+str(len(sequenceDict))
    outString += "\nTotal size: "+str(total)+"\n"

if args.data_type == "a":
    totalUngapped = sum([sum([a!="-" for a in sequenceDict[seq].seq]) for seq in sequenceDict])
    if not args.quiet_output:
        outString += "Ungapped size: "+str(totalUngapped)+" ("+str(round(totalUngapped/total*100,2))+"%)\n"
if args.calculate_base_freqs:
    outString += "Base frequencies:\n\tBase\tTotal\t(% Total"+(" / % NoGaps" if args.data_type=="a" else "")+")\n" + "\n".join(["\t"+n+"\t"+str(baseSums[n])+"\t("+str(round(baseSums[n]/total*100,2))+(("% / "+str(round(baseSums[n]/totalUngapped*100,2))) if args.data_type=="a" else "")+"%)" for n in baseSums])+"\n"
if args.calculate_n50:
    outString += "N50 sequence length: "+str(n50([len(sequenceDict[a]) for a in sequenceDict]))+"\n"
if args.calculate_pairwise_identity and args.data_type=="a":
    outString += "Mean % pairwise identity: "+str(round(sum(mpsi)/len(mpsi)*100,2))+"%\n"
elif  args.calculate_pairwise_identity and args.data_type!="a":
    print("Mean pair-wise identity can only be calculated for alignment files.")
print(outString.rstrip())