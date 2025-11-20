import pandas as pd
import os
import sys
import warnings
from helper_scripts.align_stats import *
warnings.filterwarnings("ignore")

if len(sys.argv) != 3:
    raise ValueError("Usage: python "+sys.argv[0]+" <working directory> <path to callset>")

workingdir = sys.argv[1]
path_to_call = sys.argv[2]
species = ['B.Orangutan','S.Orangutan','Gorilla','Human','Chimp','Bonobo']
allowed_chr = ['mPanTro3#1',
               'mPanPan1#P#chr1_pat_hsa1', 'mPanPan1#M#chr2_mat_hsa3', 'mPanPan1#P#chr3_pat_hsa4', 'mPanPan1#M#chr4_mat_hsa5', 'mPanPan1#P#chr5_pat_hsa6', 'mPanPan1#M#chr6_mat_hsa7', 'mPanPan1#P#chr7_pat_hsa8', 'mPanPan1#M#chr8_mat_hsa10', 'mPanPan1#P#chr9_pat_hsa11', 'mPanPan1#M#chr10_mat_hsa12', 'mPanPan1#M#chr11_mat_hsa9', 'mPanPan1#M#chr12_mat_hsa2a', 'mPanPan1#M#chr13_mat_hsa2b', 'mPanPan1#P#chr14_pat_hsa13', 'mPanPan1#M#chr15_mat_hsa14', 'mPanPan1#P#chr16_pat_hsa15', 'mPanPan1#M#chr17_mat_hsa18', 'mPanPan1#P#chr18_pat_hsa16', 'mPanPan1#P#chr19_pat_hsa17', 'mPanPan1#P#chr20_pat_hsa19', 'mPanPan1#M#chr21_mat_hsa20', 'mPanPan1#P#chr22_pat_hsa21', 'mPanPan1#P#chr23_pat_hsa22',
               'mGorGor1#P#chr1_pat_hsa1', 'mGorGor1#P#chr2_pat_hsa3', 'mGorGor1#P#chr3_pat_hsa4', 'mGorGor1#P#chr4_pat_hsa17x5', 'mGorGor1#M#chr5_mat_hsa6', 'mGorGor1#M#chr6_mat_hsa7', 'mGorGor1#P#chr7_pat_hsa8', 'mGorGor1#P#chr8_pat_hsa10', 'mGorGor1#P#chr9_pat_hsa11', 'mGorGor1#M#chr10_mat_hsa12', 'mGorGor1#M#chr11_mat_hsa2b', 'mGorGor1#P#chr12_pat_hsa2a', 'mGorGor1#P#chr13_pat_hsa9', 'mGorGor1#P#chr14_pat_hsa13', 'mGorGor1#P#chr15_pat_hsa14', 'mGorGor1#P#chr16_pat_hsa15', 'mGorGor1#M#chr17_mat_hsa18', 'mGorGor1#P#chr18_pat_hsa16', 'mGorGor1#P#chr19_pat_hsa5x17', 'mGorGor1#M#chr20_mat_hsa19', 'mGorGor1#P#chr21_pat_hsa20', 'mGorGor1#M#chr22_mat_hsa21', 'mGorGor1#M#chr23_mat_hsa22',
               'mPonPyg2#1',
               'mPonAbe1#1',
               'chm13#1']
spec_map = {'mPanTro3':'Chimp',
            'mPanPan1':'Bonobo',
            'mGorGor1':'Gorilla',
            'mPonPyg2':'B.Orangutan',
            'mPonAbe1':'S.Orangutan',
            'chm13':'Human'}

# for te in ['L1PA2','L1PA3','L1PA4']:
#     df = pd.read_csv(path_to_call+'/'+te+'.csv')
#     # filtered = df[(df['lifted_all']==1) & (df['segdup?'] == 0)]
#     for n_chr in range(1,23,1):
#         maf = dict()
#         with open(workingdir+'/chm13#1/chm13#1#chr'+str(n_chr)+'.filtered10Mb.maf','r') as c_maf:
#             for line in c_maf:
#                 if line[0] == 's':
#                     entries = line.split()
#                     chr_label_list = entries[1].split('#')
#                     # print(chr_label_list)
#                     chr_label = chr_label_list[0] + '#' + chr_label_list[1]
#                     # print(entries[1])
#                     # print(chr_label)
#                     if entries[1] in allowed_chr or chr_label in allowed_chr:
#                         # print(spec_map[entries[1].split('#')[0]])
#                         maf[entries[1]] = entries[6]
#         # pull out all annotations on this chromosome

#         # chr_filtered = filtered[filtered['query_sequence'] == 'chr'+str(n_chr)]
#         chr_filtered = df[df['query_sequence'] == 'chr'+str(n_chr)]
#         for index, row in chr_filtered.iterrows():
#             label = row['TE_label']
#             begin, end = row['query_begin']-501, row['query_end'] + 500
#             with open(workingdir+'/alignments/align_'+label+'.fasta', 'w') as fasta:
#                 # print(label)
#                 for s in maf:
#                     if s[:5] == 'chm13':
#                         fasta.write('>'+s+' chm13#1#chr'+str(n_chr)+':'+str(begin+1)+'-'+str(end)
#                                     +'\n'+maf[s][begin:end]+'\n')
#                     else:
#                         # print(f'{label} {s} {n_chr}')
#                         # print()
#                         seq = maf[s][begin:end]
#                         fasta.write('>'+s+'\n'+seq+'\n') 


alignment_stats_wfmash(workingdir)




