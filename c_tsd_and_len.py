import sys
import pandas as pd
from helper_scripts.pull_seq import pull
from helper_scripts.TSD import * 

if len(sys.argv) != 2:
    raise ValueError("Usage: python "+sys.argv[0]+" <working directory>")

path = sys.argv[1] # working directory path
genome = pull(path+'/genomes/hs1.fa')

# pull out the sequnce to check for tsd in
te_list = ['L1PA2','L1PA3','L1PA4']
num_tsd = [0,0,0]
num_full_len = [0,0,0]
for te in te_list:
    df = pd.read_csv(path+'/call_set/'+te+'.csv')
    df['TSD_found'], df['full_len'] = 0, 0
    df['TSD_seq'], df['TSD_start'], df['TSD_end'] = None, None, None
    for index, row in df.iterrows():
        start1, start2 = row['query_begin']-51, row['query_end']-15
        tsd_seq1 = genome[row['query_sequence']][row['query_begin']-51:row['query_begin']+14]
        tsd_seq2 = genome[row['query_sequence']][row['query_end']-15:row['query_end']+50] 
        seq, off1, off2 = tsd_human(tsd_seq1,tsd_seq2)
        if seq:
            df.loc[index, 'TSD_found'] = 1
            df.loc[index, 'TSD_seq'] = seq
            df.loc[index, 'TSD_start'] = start1 + off1
            df.loc[index, 'TSD_end'] = start2 + off2
            num_tsd[te_list.index(te)] += 1
        if row['region_length'] > 5900:
            df.loc[index,'full_len'] = 1
            num_full_len[te_list.index(te)] += 1
    df.to_csv(path+'/call_set/'+te+'.csv', index=False)
    
print(f'number of tsds found in human {num_tsd[0]}, {num_tsd[1]}, {num_tsd[2]}')
print(f'number of full length elements {num_full_len[0]}, {num_full_len[1]}, {num_full_len[2]}')


