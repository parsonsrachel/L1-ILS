
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

orf_df = read_orf_file(workingdir)
print(orf_df['L1PA2']['orf1_ref'].loc[orf_df['L1PA2']['orf1_ref']['TE label'] == 'L1PA2_0'])

columns=['B.Orangutan','S.Orangutan','Gorilla','Human','Chimp','Bonobo']

# add data to the call_set df
for k in dfs:
    temp_df = pd.read_csv(workingdir+'/call_set/'+k+'.csv')
    # for c in columns:
    #     temp_df = temp_df.drop(c+'_seq',axis=1)
    temp_df[['orthologous','presence/absence','TSD_pattern','ORF1','ORF2','gene_tree']] = [None,None,None,None,None,None]
    
    tsd_df = pd.DataFrame(columns=['TE_label','species','presence/absence','TSD1','TSD2','TSD_pattern'])
    for index, row in temp_df.iterrows():
        label = row['TE_label']
        pattern = dict()
        orth = dict()
        tsd_pattern = dict()
        orf1,orf2 = None,None
        for c in columns:
            # print(dfs[k].loc[(dfs[k]['TE label'] == label) and (dfs[k]['species'] == c),['orthologous','presence/absence']])
            # o, pa = dfs[k][(dfs[k]['TE label'] == label) & (dfs[k]['species'] == c),['orthologous','presence/absence']].iloc[0]
            if c == 'Human':
                orth[c] = '1'
                pattern[c] = '1'
                
            else:
                filtered = dfs[k][(dfs[k]['TE label'] == label) & (dfs[k]['species'] == c)]
                if not filtered.empty:
                    o, pa = filtered[['orthologous', 'presence/absence']].iloc[0]
                else:
                    o, pa = None, None  # or handle it however you'd like

                orth[c] = o
                pattern[c] = pa
            tsd_pattern[c] = None
        
        if row['TSD_found'] and row['lifted_all'] and not row['segdup?']:
            begin,seq, tsd1_start, tsd2_start = row[['query_begin','TSD_seq','TSD_start','TSD_end']]
            msa = read_align(workingdir+'/alignments/align_'+label+'.fasta')
            # print(msa)

            res = tsd_non_human(msa,(tsd1_start - begin) + 501,(tsd2_start - begin) + 501,seq)
            # print(res)
            # get the orth, pattern, tsd pattern, orf1, orf2 calls
            for c in columns:
                if c == 'Human':
                    tsd_pattern[c] = '1'
                else:
                    tsd_pattern[c] = res[c][2]
                tsd_df.loc[len(tsd_df)] = [label,c,pattern[c],res[c][0],res[c][1],res[c][2]]
        
        gene_tree = ''
        try:
            with open(workingdir+'/trees/te/align_'+label+'.fasta.treefile','r') as treefile:
                for line in treefile:
                    gene_tree+= line.strip()
        except:
            gene_tree = ''
        
        w_orth, w_pattern, w_tsd = [],[],[]
        w_orf1, w_orf2 = [], []
        for c in columns:
            w_orth.append(orth[c])
            w_pattern.append(pattern[c])
            w_tsd.append(tsd_pattern[c])
            # print(label)
            orf1_row = orf_df[k]['orf1_ref'].loc[(orf_df[k]['orf1_ref']['TE label'] == label) & (orf_df[k]['orf1_ref']['species'] == c)]
            orf2_row = orf_df[k]['orf2_ref'].loc[(orf_df[k]['orf2_ref']['TE label'] == label) & (orf_df[k]['orf2_ref']['species'] == c)]
            # print(orf1_row.shape[0])
            # print(type(orf2_row.shape[0]))
            if orf1_row.shape[0] != 0:
                for _, row in orf1_row.iterrows():
                    orf1 = orf_stats('orf1_ref',row)
            else:
                orf1 = [0]
            
            if orf2_row.shape[0] != 0:
                for _, row in orf2_row.iterrows():
                    orf2 = orf_stats('orf2_ref',row)
            else:
                orf2 = [0]
            w_orf1.append(orf1[0])
            w_orf2.append(orf2[0])
        # print(w_orth)
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

            str_orf1 += comma + str(w_orf1[i])
            str_orf2 += comma + str(w_orf2[i])


        # print(w_pattern)
        # print(w_tsd)
        print(f'label:{label} gene tree: {gene_tree}')
        temp_df.loc[index,['orthologous','presence/absence','TSD_pattern','ORF1','ORF2','gene_tree']] =[str_orth,str_pat,str_tsd,str_orf1,str_orf2,gene_tree]
    temp_df.to_csv(workingdir+'/call_set/'+k+'.csv',index=False)
    tsd_df.to_csv(workingdir+'/tsd/'+k+'.csv',index=False)


