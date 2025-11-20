import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

def filter(df):
    df['fragments_length'] = df['query_end'] - df['query_begin'] + 1
    # print(df.sort_values(by=['query_sequence','query_begin']).head(10))

    df = df.groupby(['repeat_name', 'ID']).agg({
        "SW_score": lambda x: '/'.join(str(i) for i in x), 
        "perc_div": lambda x: '/'.join(str(i) for i in x), 
        "perc_del": lambda x: '/'.join(str(i) for i in x), 
        "perc_ins": lambda x: '/'.join(str(i) for i in x),
        'query_sequence': 'first',       
        'query_begin': 'first',
        'query_end': 'last',
        "query_left": lambda x: '/'.join(x),
        'strand': lambda x: '/'.join(x),   
        "repeat_class_family": 'first',
        "repeat_begin": lambda x: '/'.join(x), 
        "repeat_end": lambda x: '/'.join(str(i) for i in x), 
        "repeat_left": lambda x: '/'.join(x),
        "fragments_length": lambda x: sum(x)
    }).reset_index()

    # add in the length of the TE region and the total # of bp between aligned regions (some are negative because there are overlapping annotations)
    df['region_length'] = df['query_end'] - df['query_begin'] + 1
    df['dist_between_regions'] = df['region_length'] - df['fragments_length']

    # get counts for each element
    total_counts = [0,0,0]
    te_list = ['L1PA2','L1PA3','L1PA4']
    for index, row in df.iterrows():
        total_counts[te_list.index(row['repeat_name'])]+=1

    # print(f'In total there are {total_counts} of L1PA2,3,4')

    # filter out repeats with region_length < 200 
    df = df[df['region_length'] >= 200]

    first_filt_counts = [0,0,0]
    for index, row in df.iterrows():
        first_filt_counts[te_list.index(row['repeat_name'])]+=1

    # print(f'In total there are {first_filt_counts} of L1PA2,3,4 after first filter (len >= 200)')
    
    # filter out repeats less than 500 bp to another
    temp = df.sort_values(by=['query_sequence','query_begin']).reset_index(drop=True)
    keep = [True] * len(temp)

    for i in range(len(temp) - 1):
        curr = temp.loc[i, 'query_end']
        future = temp.loc[i + 1, 'query_begin']
        if (future - curr) < 500:
            keep[i] = False
            keep[i+1] = False
    df = temp[keep].reset_index(drop=True)

    sec_filt_counts = [0,0,0]
    for index, row in df.iterrows():
        sec_filt_counts[te_list.index(row['repeat_name'])]+=1

    # print(f'In total there are {sec_filt_counts} of L1PA2,3,4')
    
    return df, total_counts, first_filt_counts, sec_filt_counts
