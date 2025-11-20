
import pandas as pd
import sys
from helper_scripts.align_stats import *
from helper_scripts.TSD import *
from helper_scripts.orf_stats import *

if len(sys.argv) != 2:
    raise ValueError("Usage: python "+sys.argv[0]+" <working directory>")

workingdir = sys.argv[1] # working directory path

# dfs = alignment_stats(workingdir)
dfs = dict()
for i in ['2','3','4']:
    dfs['L1PA'+i] = pd.read_csv(workingdir+'/alignment_data/output'+i+'.csv')


columns=['B.Orangutan','S.Orangutan','Gorilla','Human','Chimp','Bonobo']
prefixes = {'mPonP':'B.Orangutan',
            'mPonA':'S.Orangutan',
            'mGorG':'Gorilla',
            'Human':'Human',
            'mPanT':'Chimp',
            'mPanP':'Bonobo'}

# add data to the call_set df
for k in dfs:
    temp_df = pd.read_csv(workingdir+'/call_set/'+k+'.csv')
    # for c in columns:
    #     temp_df = temp_df.drop(c+'_seq',axis=1)
    temp_df[['wf_orthologous','wf_presence/absence','wf_TSD_pattern']] = [None,None,None]
    
    tsd_df = pd.DataFrame(columns=['TE_label','species','presence/absence','TSD1','TSD2','TSD_pattern'])
    for index, row in temp_df.iterrows():
        label = row['TE_label']
        pattern = dict()
        orth = dict()
        tsd_pattern = dict()
        filtered = dfs[k][(dfs[k]['TE label'] == label)]
        print(filtered.head())
        msa_labels = []
        
        if filtered.shape[0] == 6:
            for d_index, d_row in filtered.iterrows():
                msa_labels.append(d_row['species'])
                c = d_row['species'][:5]
                orth[prefixes[c]] = d_row['orthologous']
                pattern[prefixes[c]] = d_row['presence/absence']
                tsd_pattern[prefixes[c]] = None
        else:
            for p in prefixes:
                dup_filt = filtered[filtered['species'].str[:5] == p]
                if dup_filt.shape[0] > 1:
                    # take whichever establishes orthology
                    good_ind = 0
                    for d_index, d_row in dup_filt.iterrows():
                        if d_row['orthologous'] == 'TRUE':
                            good_ind = d_index
                            break
                    c = dup_filt.loc[d_index,'species'][:5]
                    msa_labels.append(dup_filt.loc[d_index,'species'])
                    orth[prefixes[c]] = dup_filt.loc[d_index,'orthologous']
                    pattern[prefixes[c]] = dup_filt.loc[d_index,'presence/absence']
                    tsd_pattern[prefixes[c]] = None
                else:
                    c = dup_filt['species'].iloc[0][:5]
                    msa_labels.append(dup_filt['species'].iloc[0])
                    orth[prefixes[c]] = dup_filt['orthologous'].iloc[0]
                    pattern[prefixes[c]] = dup_filt['presence/absence'].iloc[0]
                    tsd_pattern[prefixes[c]] = None 
            
        # print('orth dictionary:')
        # print(orth)
        if row['TSD_found']:
            begin,seq, tsd1_start, tsd2_start = row[['query_begin','TSD_seq','TSD_start','TSD_end']]
            msa = read_align(workingdir+'/alignments/align_'+label+'.fasta')
            print(msa.keys())
            for key in msa:
                if key[:5] == 'chm13':
                    msa['Human'] = msa.pop(key)
                    break
            

            res = tsd_non_human(msa,(tsd1_start - begin) + 501,(tsd2_start - begin) + 501,seq)
            print(res)
            print(msa_labels)
            # get the orth, pattern, tsd pattern, orf1, orf2 calls
            for key in msa_labels:
                if key == 'Human':
                    tsd_pattern['Human'] = '1'
                else:
                    c = prefixes[key[:5]]
                    tsd_pattern[c] = res[key][2]
                    tsd_df.loc[len(tsd_df)] = [label,c,pattern[c],res[key][0],res[key][1],res[key][2]]
        
        w_orth, w_pattern, w_tsd = [],[],[]
        w_orf1, w_orf2 = [], []
        for c in columns:
            w_orth.append(orth[c])
            w_pattern.append(pattern[c])
            w_tsd.append(tsd_pattern[c])

        str_orth = ''
        str_pat = ''
        str_tsd = ''
        str_orf1 = ''
        str_orf2 = ''
        comma = ''
        for i in range(len(w_orth)):
            if i != 0:
                comma = ','

            if w_orth[i] is None:
                str_orth += comma + 'None' 
            if w_orth[i] == '1' or w_orth[i] == 1:
                str_orth += comma + '1'
            elif w_orth[i] == '0' or w_orth[i] == 0:
                str_orth += comma + '0'

            if w_pattern[i] is None:
                str_pat += comma + 'None' 
            if w_pattern[i] == '1' or w_pattern[i] == 1:
                str_pat += comma + '1'
            elif w_pattern[i] == '0' or w_pattern[i] == 0:
                str_pat += comma + '0'
            elif w_pattern[i] == '?':
                str_pat += comma + '?' 

            if w_tsd[i] is None:
                str_tsd += comma + 'None' 
            if w_tsd[i] == '1' or w_tsd[i] == 1:
                str_tsd += comma + '1'
            elif w_tsd[i] == '0' or w_tsd[i] == 0:
                str_tsd += comma + '0'
            elif w_tsd[i] == '?':
                str_tsd += comma + '?' 

            # str_orf1 += comma + str(w_orf1[i])
            # str_orf2 += comma + str(w_orf2[i])

        print(label)
        print(str_pat)
        print(str_orth)

        temp_df.loc[index,['wf_orthologous','wf_presence/absence','wf_TSD_pattern']] =[str_orth,str_pat,str_tsd]
    temp_df.to_csv(workingdir+'/call_set/'+k+'.csv',index=False)
    tsd_df.to_csv(workingdir+'/tsd/'+k+'.csv',index=False)


