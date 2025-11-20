import pandas as pd
import os
import itertools
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

class msaStats:
    def __init__(self, msa, te_regions, begin):
        self.data = dict()
        self.te_regions = te_regions
        self.begin = begin
        self.msa = msa

        # check that the file was read correctly 
        # and each sequence is the same length
        labels = msa.keys()
        for i in itertools.combinations(labels, 2):
            assert(len(msa[i[0]]) == len(msa[i[1]]))

        apes = ['Chimp','S.Orangutan','B.Orangutan','Bonobo','Gorilla']
        # apes = ['Chimp','Bonobo']
        for a in apes:
            self.data[a] = self._spec_stats(a)

    def both_gaps(self,a):
        return self.data[a][0]
    
    def gap_in_human(self,a):
        return self.data[a][1]
    
    def gap_in_ape(self,a):
        return self.data[a][2]
    
    def no_gaps(self,a):
        return self.data[a][3]
    
    def nuc_match(self,a):
        return self.data[a][4]
    
    def _spec_stats(self,a):
        both_gap = [0, 0, 0]
        gap_h = [0, 0, 0]
        gap_a = [0, 0, 0]
        no_gap = [0, 0, 0]
        same_nuc = [0, 0, 0]
        seq_ind = -1
        prev_ind = -1
        len_human = len(self.msa['Human'].replace('-',''))
        # print(msa['Human'])
        for i in range(len(self.msa['Human'])):
            if self.msa['Human'][i] != '-':
                prev_ind = seq_ind
                seq_ind += 1

            if self.msa['Human'][i] == '-' and self.msa[a][i] == '-':
                #gap in both
                if seq_ind < 500:
                    both_gap[0] += 1
                elif seq_ind >= len_human - 500:
                    both_gap[1] += 1
                elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                    both_gap[2] += 1

            elif self.msa['Human'][i] == '-' and self.msa[a][i] != '-':
                # gap only in human
                if seq_ind < 500 and prev_ind == seq_ind:
                    gap_h[0] += 1
                elif seq_ind >= len_human - 500:
                    gap_h[1] += 1
                elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                    gap_h[2] += 1

            elif self.msa['Human'][i] != '-' and self.msa[a][i] == '-':
                # gap only in other ape
                if seq_ind < 500:
                    gap_a[0] += 1
                elif seq_ind >= len_human - 500:
                    gap_a[1] += 1
                elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                    gap_a[2] += 1
                    
            elif self.msa['Human'][i] != '-' and self.msa[a][i] != '-':
                #gap in neither
                if self.msa['Human'][i].upper() == self.msa[a][i].upper():
                    # same nucleotide in both
                    if seq_ind < 500:
                        same_nuc[0] += 1
                    elif seq_ind >= len_human - 500:
                        same_nuc[1] += 1
                    elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                        same_nuc[2] += 1
                if seq_ind < 500:
                    no_gap[0] += 1
                elif seq_ind >= len_human - 500:
                    no_gap[1] += 1
                elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                    no_gap[2] += 1

        return [both_gap,gap_h,gap_a,no_gap,same_nuc]

class msaStatsWF:
    def __init__(self, msa, te_regions, begin, chroms):
        self.data = dict()
        self.te_regions = te_regions
        self.begin = begin
        self.msa = msa

        # check that the file was read correctly 
        # and each sequence is the same length
        labels = msa.keys()
        for i in itertools.combinations(labels, 2):
            assert(len(msa[i[0]]) == len(msa[i[1]]))

        # apes = ['Chimp','S.Orangutan','B.Orangutan','Bonobo','Gorilla']
        # apes = ['Chimp','Bonobo']
        for a in chroms:
            self.data[a] = self._spec_stats(a)

    def both_gaps(self,a):
        return self.data[a][0]
    
    def gap_in_human(self,a):
        return self.data[a][1]
    
    def gap_in_ape(self,a):
        return self.data[a][2]
    
    def no_gaps(self,a):
        return self.data[a][3]
    
    def nuc_match(self,a):
        return self.data[a][4]
    
    def _spec_stats(self,a):
        both_gap = [0, 0, 0]
        gap_h = [0, 0, 0]
        gap_a = [0, 0, 0]
        no_gap = [0, 0, 0]
        same_nuc = [0, 0, 0]
        seq_ind = -1
        prev_ind = -1
        len_human = len(self.msa['Human'].replace('-',''))
        # print(msa['Human'])
        for i in range(len(self.msa['Human'])):
            if self.msa['Human'][i] != '-':
                prev_ind = seq_ind
                seq_ind += 1

            if self.msa['Human'][i] == '-' and self.msa[a][i] == '-':
                #gap in both
                if seq_ind < 500:
                    both_gap[0] += 1
                elif seq_ind >= len_human - 500:
                    both_gap[1] += 1
                elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                    both_gap[2] += 1

            elif self.msa['Human'][i] == '-' and self.msa[a][i] != '-':
                # gap only in human
                if seq_ind < 500 and prev_ind == seq_ind:
                    gap_h[0] += 1
                elif seq_ind >= len_human - 500:
                    gap_h[1] += 1
                elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                    gap_h[2] += 1

            elif self.msa['Human'][i] != '-' and self.msa[a][i] == '-':
                # gap only in other ape
                if seq_ind < 500:
                    gap_a[0] += 1
                elif seq_ind >= len_human - 500:
                    gap_a[1] += 1
                elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                    gap_a[2] += 1
                    
            elif self.msa['Human'][i] != '-' and self.msa[a][i] != '-':
                #gap in neither
                if self.msa['Human'][i].upper() == self.msa[a][i].upper():
                    # same nucleotide in both
                    if seq_ind < 500:
                        same_nuc[0] += 1
                    elif seq_ind >= len_human - 500:
                        same_nuc[1] += 1
                    elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                        same_nuc[2] += 1
                if seq_ind < 500:
                    no_gap[0] += 1
                elif seq_ind >= len_human - 500:
                    no_gap[1] += 1
                elif ((self.te_regions['query_begin'] <= self.begin + seq_ind) & (self.te_regions['query_end'] >= self.begin + seq_ind)).any():
                    no_gap[2] += 1

        return [both_gap,gap_h,gap_a,no_gap,same_nuc]
    
def alignment_stats(path):
    col_labels = ['TE label','species','length of TE region in human', 'orthologous', 'presence/absence',
                'sequence similarity f1', 'sequence similarity f2', 'sequence similarity TE',
                'both gaps f1','both gaps f2','both gaps TE',
                'gap in human not ape f1','gap in human not ape f2','gap in human not ape TE',
                'gap in ape not human f1','gap in ape not human f2','gap in ape not human TE',
                'both nongap f1','both nongap f2','both nongap TE',
                'exact match f1','exact match f2','exact match TE'
                ]


    rm_df = pd.read_csv(path+'/repeat_masker_out/rm_annotations.out',
                    sep=r'\s+', 
                    names=[
                        "SW_score", "perc_div", "perc_del", "perc_ins",
                        "query_sequence", "query_begin", "query_end", "query_left",
                        "strand", "repeat_name", "repeat_class_family",
                        "repeat_begin", "repeat_end", "repeat_left", "ID","star"
                    ])
    output_dfs = dict()
    for l in ['2','3','4']:
        df = pd.DataFrame(columns = col_labels)
        for filename in os.listdir(path+'/alignments/'):
            if filename != 'intemp.fa' and filename.split('_')[1][-1] == l:
                msa = dict()
                label = ''
                chr= ''
                begin, end = 0, 0

                with open(path+'/alignments/'+filename,'r') as fasta_file:
                    spec = ''
                    for line in fasta_file:
                        if line[0] == '>':
                            # save the species label
                            sp = line.split()
                            spec = sp[0][1:]
                            if spec == 'Human':
                                chr = sp[1].split('_')[0].split(':')[0].split('#')[-1]
                                begin, end = sp[1].split('_')[0].split(':')[1].split('-')
                                begin, end = int(begin), int(end)
                                # print(chr)
                                # print(begin, end)
                            if label != '':
                                assert(label == sp[1].split('_')[-3] + '_' + sp[1].split('_')[-2])
                            else:
                                label = sp[1].split('_')[-3] + '_' + sp[1].split('_')[-2]
                            
                            # print(label)
                            msa[spec] = ''
                        else:
                            msa[spec] += line.strip()

                id = rm_df.loc[(rm_df['query_sequence'] == chr) & (rm_df['query_begin'] == begin + 500), 'ID']
                te_regions = rm_df[rm_df['ID'].isin(id)][['query_begin','query_end']]
                len_human = 0

                for _, row in te_regions.iterrows():
                    len_human += row['query_end'] - row['query_begin'] + 1

                if len_human == 0:
                    print(filename)

                stats = msaStats(msa, te_regions, begin)

                for a in ['Chimp','S.Orangutan','B.Orangutan','Bonobo','Gorilla']:
                    seq_sim = [(x / y if x != 0 else 0) for x, y in zip(stats.nuc_match(a),stats.no_gaps(a))]

                    orth = False
                    no_gaps = stats.no_gaps(a)
                    if no_gaps[0] > 250 and no_gaps[1] > 250 and seq_sim[0] > 0.9 and seq_sim[1] > 0.9:
                        # both flanking regions are less than 50% gaps
                        orth = True

                    if no_gaps[2]/len_human < 0.1:
                        pattern = '0'
                    elif no_gaps[2]/len_human > 0.9 and seq_sim[2] > 0.9:
                        pattern = '1'
                    else:
                        pattern = '?'

                    list1 = [label, a, len_human, orth, pattern]

                    combined_iterator = itertools.chain(list1, 
                                                        seq_sim,
                                                        stats.both_gaps(a), 
                                                        stats.gap_in_human(a), 
                                                        stats.gap_in_ape(a), 
                                                        stats.no_gaps(a), 
                                                        stats.nuc_match(a))
                    combined_list = list(combined_iterator)
                    df.loc[len(df)] = combined_list

        df.to_csv(path+"/alignment_data/output"+l+".csv", index=False)
        output_dfs['L1PA'+l] = df
        del df
    
    return output_dfs

def alignment_stats_wfmash(path):
    col_labels = ['TE label','species','length of TE region in human', 'orthologous', 'presence/absence',
                'sequence similarity f1', 'sequence similarity f2', 'sequence similarity TE',
                'both gaps f1','both gaps f2','both gaps TE',
                'gap in human not ape f1','gap in human not ape f2','gap in human not ape TE',
                'gap in ape not human f1','gap in ape not human f2','gap in ape not human TE',
                'both nongap f1','both nongap f2','both nongap TE',
                'exact match f1','exact match f2','exact match TE'
                ]


    rm_df = pd.read_csv(path+'/repeat_masker_out/rm_annotations.out',
                    sep=r'\s+', 
                    names=[
                        "SW_score", "perc_div", "perc_del", "perc_ins",
                        "query_sequence", "query_begin", "query_end", "query_left",
                        "strand", "repeat_name", "repeat_class_family",
                        "repeat_begin", "repeat_end", "repeat_left", "ID","star"
                    ])
    output_dfs = dict()
    for l in ['2','3','4']:
        df = pd.DataFrame(columns = col_labels)
        for filename in os.listdir(path+'/alignments/'):
            if filename != 'intemp.fa' and filename.split('_')[1][-1] == l:
                msa = dict()
                label = ''
                chr= ''
                begin, end = 0, 0

                with open(path+'/alignments/'+filename,'r') as fasta_file:
                    spec = ''
                    label = filename.split('.')[0].split('_',1)[1]
                    for line in fasta_file:
                        if line.startswith('>'):
                            sp = line.split()
                            spec = sp[0][1:]
                            if line[1:6] == 'chm13':
                                chr = sp[1].split('#')[-1].split(':')[0]
                                begin, _ = map(int, sp[1].split(':')[-1].split('-'))
                                label = filename.split('.')[0].split('_',1)[1]
                                msa['Human'] = ''
                                spec = 'Human'
                            else:
                                msa[spec] = ''
                        else:
                            msa[spec] += line.strip()
                print(msa.keys())
                id = rm_df.loc[(rm_df['query_sequence'] == chr) & (rm_df['query_begin'] == begin + 500), 'ID']
                te_regions = rm_df[rm_df['ID'].isin(id)][['query_begin','query_end']]
                len_human = 0
                # print(te_regions.head())
                # print(msa['Bonobo'][:500])
                # print(msa['Human'][:500])

                for _, row in te_regions.iterrows():
                    len_human += row['query_end'] - row['query_begin'] + 1

                if len_human == 0:
                    print(f'{chr} {begin} {end}')
                    print(filename)

                stats = msaStatsWF(msa, te_regions, begin, list(msa.keys()))

                for a in msa:
                    seq_sim = [(x / y if x != 0 else 0) for x, y in zip(stats.nuc_match(a),stats.no_gaps(a))]

                    orth = False
                    no_gaps = stats.no_gaps(a)
                    if no_gaps[0] > 250 and no_gaps[1] > 250 and seq_sim[0] > 0.9 and seq_sim[1] > 0.9:
                        # both flanking regions are less than 50% gaps
                        orth = True

                    if no_gaps[2]/len_human < 0.1:
                        pattern = '0'
                    elif no_gaps[2]/len_human > 0.9 and seq_sim[2] > 0.9:
                        pattern = '1'
                    else:
                        pattern = '?'

                    list1 = [label, a, len_human, orth, pattern]

                    combined_iterator = itertools.chain(list1, 
                                                        seq_sim,
                                                        stats.both_gaps(a), 
                                                        stats.gap_in_human(a), 
                                                        stats.gap_in_ape(a), 
                                                        stats.no_gaps(a), 
                                                        stats.nuc_match(a))
                    combined_list = list(combined_iterator)
                    df.loc[len(df)] = combined_list

        df.to_csv(path+"/alignment_data/output"+l+".csv", index=False)
        output_dfs['L1PA'+l] = df
        del df
        print(f'finished L1PA{l}')
    
    return output_dfs

# def process_alignment_file(args):
#     filename, l, path, rm_df = args

#     col_labels = ['TE label','species','length of TE region in human', 'orthologous', 'presence/absence',
#                   'sequence similarity f1', 'sequence similarity f2', 'sequence similarity TE',
#                   'both gaps f1','both gaps f2','both gaps TE',
#                   'gap in human not ape f1','gap in human not ape f2','gap in human not ape TE',
#                   'gap in ape not human f1','gap in ape not human f2','gap in ape not human TE',
#                   'both nongap f1','both nongap f2','both nongap TE',
#                   'exact match f1','exact match f2','exact match TE']

#     df = pd.DataFrame(columns = col_labels)

#     msa = dict()
#     label = ''
#     chr = ''
#     begin, end = 0, 0

#     with open(path+'/alignments/'+filename, 'r') as fasta_file:
#         for line in fasta_file:
#             if line.startswith('>'):
#                 sp = line.split()
#                 spec = sp[0][1:]
#                 if line[:6] == 'chm13':
#                     chr = sp[1].split('#')[-1].split(':')[0]
#                     begin, _ = map(int, sp[1].split(':')[-1].split('-'))
#                     label = filename.split('.')[0].split('_',1)[1]
#                 msa[spec] = ''
#             else:
#                 msa[spec] += line.strip()

#     id = rm_df.loc[(rm_df['query_sequence'] == chr) & (rm_df['query_begin'] == begin + 500), 'ID']
#     te_regions = rm_df[rm_df['ID'].isin(id)][['query_begin','query_end']]
#     len_human = sum(row['query_end'] - row['query_begin'] + 1 for _, row in te_regions.iterrows())

#     if len_human == 0:
#         return df  # skip empty results

#     stats = msaStats(msa, te_regions, begin)

#     for a in msa:
#         seq_sim = [(x / y if y != 0 else 0) for x, y in zip(stats.nuc_match(a),stats.no_gaps(a))]

#         orth = False
#         no_gaps = stats.no_gaps(a)
#         if no_gaps[0] > 250 and no_gaps[1] > 250 and seq_sim[0] > 0.9 and seq_sim[1] > 0.9:
#             orth = True

#         if no_gaps[2]/len_human < 0.1:
#             pattern = '0'
#         elif no_gaps[2]/len_human > 0.9 and seq_sim[2] > 0.9:
#             pattern = '1'
#         else:
#             pattern = '?'

#         list1 = [label, a, len_human, orth, pattern]
#         combined_list = list(itertools.chain(list1, seq_sim,
#                                              stats.both_gaps(a), 
#                                              stats.gap_in_human(a), 
#                                              stats.gap_in_ape(a), 
#                                              stats.no_gaps(a), 
#                                              stats.nuc_match(a)))
#         df.loc[len(df)] = combined_list

#     return df

# def alignment_stats_wfmash(path):
#     col_labels = ['TE label','species','length of TE region in human', 'orthologous', 'presence/absence',
#                   'sequence similarity f1', 'sequence similarity f2', 'sequence similarity TE',
#                   'both gaps f1','both gaps f2','both gaps TE',
#                   'gap in human not ape f1','gap in human not ape f2','gap in human not ape TE',
#                   'gap in ape not human f1','gap in ape not human f2','gap in ape not human TE',
#                   'both nongap f1','both nongap f2','both nongap TE',
#                   'exact match f1','exact match f2','exact match TE']

#     # rm_df = pd.read_csv(path+'/repeat_masker_out/rm_annotations.out',
#     #     sep=r'\s+', 
#     #     names=[...])  # same as before
#     rm_df = pd.read_csv(path+'/repeat_masker_out/rm_annotations.out',
#                     sep=r'\s+', 
#                     names=[
#                         "SW_score", "perc_div", "perc_del", "perc_ins",
#                         "query_sequence", "query_begin", "query_end", "query_left",
#                         "strand", "repeat_name", "repeat_class_family",
#                         "repeat_begin", "repeat_end", "repeat_left", "ID","star"
#                     ])

#     output_dfs = dict()

#     for l in ['2','3','4']:
#         # Get filenames that match the current group
#         filenames = [f for f in os.listdir(path+'/alignments/') 
#                      if f != 'intemp.fa' and f.split('_')[1][-1] == l]

#         # Build argument list
#         args = [(filename, l, path, rm_df) for filename in filenames]

#         # Use multiprocessing to parallelize
#         with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
#             dfs = list(executor.map(process_alignment_file, args))

#         # Combine all results into one DataFrame
#         combined_df = pd.concat(dfs, ignore_index=True)
#         combined_df.to_csv(path+f"/alignment_data/output{l}.csv", index=False)
#         output_dfs['L1PA'+l] = combined_df

#     return output_dfs

def read_align(filename):
    msa = dict()
    with open(filename,'r') as fa:
        spec = ''
        for line in fa:
            if line[0] == '>':
                spec = line.split()[0][1:]
                
            else:
                if spec in msa:
                    msa[spec] += line.strip()
                else:
                    msa[spec] = line.strip()
    return msa

# alignment_stats_wfmash('/Users/rachelparsons/Downloads/T2T_primates/L1PA_ILS_Analysis/wfmash')
