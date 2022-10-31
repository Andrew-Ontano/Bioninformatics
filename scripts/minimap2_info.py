import sys
import pandas as pd


header = ['QueryName', 'QueryLength', 'QueryStart', 'QueryEnd', 'QueryStrand', 'TargetName', 'TargetLength', 'TargetStart', 'TargetEnd', 'MapBases', 'MapSize', 'MapQuality', 'AlignmentType', 'ChainMinimizers', 'Chain1Score', 'Chain2Score', 'SequenceDivergence', 'RepetitiveSeedLength']
dataframe = pd.read_csv(sys.argv[1], sep='\t', names=header)

#for group in dataframe.groupby(['QueryName']):
    #tempMatch = group[1].head(n=1)
    #print(f"{tempMatch.iloc[0, 0]}\t{tempMatch.iloc[0, 5]}")

dataframe = dataframe.loc[(dataframe['MapBases'] / dataframe['QueryLength']) > 0.5]
dataframe.to_csv(sys.argv[2], sep='\t', header=False, index=False)

