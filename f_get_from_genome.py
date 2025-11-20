import pandas as pd
import sys
import warnings
from helper_scripts.pull_seq import *
warnings.filterwarnings("ignore")


if len(sys.argv) != 2:
    raise ValueError("Usage: python "+sys.argv[0]+" <working directory>")

path = sys.argv[1] # working directory path

genome_map = {path+'/genomes/hs1.fa':'Human',
                path+'/genomes/mPonPyg2.dip.cur.20231122.fasta':'B.Orangutan',
                path+'/genomes/mPonAbe1.dip.cur.20231205.fasta':'S.Orangutan',
                path+'/genomes/mPanTro3.dip.cur.20231122.fasta':'Chimp',
                path+'/genomes/mGorGor1.dip.cur.20231122.fasta':'Gorilla',
                path+'/genomes/mPanPan1.dip.cur.20231122.fasta':'Bonobo'}
tes = ['L1PA2','L1PA3','L1PA4']
columns=['B.Orangutan','B.Orangutan seq','S.Orangutan','S.Orangutan seq','Gorilla','Gorilla seq','Human','Human seq','Chimp','Chimp seq','Bonobo','Bonobo seq']


out_csvs = {'L1PA2':pd.DataFrame(columns=columns),'L1PA3':pd.DataFrame(columns=columns),'L1PA4':pd.DataFrame(columns=columns)}
for s_filename in genome_map:
    read_genome = pull(s_filename)
    taxa =genome_map[s_filename]
    
    # loop over te family
    for te in tes:
        call_set_df = pd.read_csv(path+'/call_set/'+te+'.csv')
        call_set_df[taxa+'_seq'] = ''
        coord_file_map = {'Human':path+'/coords_to_pull/'+te+'/Human.csv',
                        'B.Orangutan':path+'/coords_to_pull/'+te+'/B.Orangutan.csv',
                        'S.Orangutan':path+'/coords_to_pull/'+te+'/S.Orangutan.csv',
                        'Chimp':path+'/coords_to_pull/'+te+'/Chimp.csv',
                        'Gorilla':path+'/coords_to_pull/'+te+'/Gorilla.csv',
                        'Bonobo':path+'/coords_to_pull/'+te+'/Bonobo.csv'}
        # loop over coordinates
        out_sequences = []
        coord_df = pd.read_csv(coord_file_map[taxa],header=None)
        coord_df.columns = ['chr:begin-end']

        for index,row in coord_df.iterrows():
            # print(row[0])
            if taxa == 'Human':
                #read file different because it has different format
                chr, begin, end, label = row[0].split()
                strand = '+'
                # chr = 'chm13#1#' + chr
            else:
                chr = row[0].split('#')[-1].split(':')[0]
                strand = row[0].strip().split('_')[-1]
                label = '_'.join(row[0].strip().split('#')[-1].split('_')[-4:][:2]) 
                begin, end = row[0].strip().split('#')[-1].split(':')[-1].split('_')[0].split('-')

            begin, end = int(begin), int(end)
            chr_str = chr + ':' + str(begin) + '-' + str(end) + '_' + label + '_' + strand
            # print(chr_str)
            if strand == '-':
                call_set_df.loc[(call_set_df['TE_label'] == label),taxa+'_seq'] = rv_comp(read_genome[chr][begin-1:end]) 
                out_csvs[te].loc[index, [taxa, taxa+' seq']] = [chr_str, rv_comp(read_genome[chr][begin-1:end])]
            else:
                call_set_df.loc[(call_set_df['TE_label'] == label),taxa+'_seq'] = read_genome[chr][begin-1:end]
                out_csvs[te].loc[index, [taxa, taxa+' seq']] = [chr_str, read_genome[chr][begin-1:end]]

        call_set_df.to_csv(path+'/call_set/'+te+'.csv',index=False)

for te in tes:
    out_csvs[te].to_csv(path+'/sequences/'+te+'.csv', sep='\t')


    
