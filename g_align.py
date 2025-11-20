# script to align the sequences fould using get_from_genome.py
import pandas as pd
import os
import sys
import warnings
warnings.filterwarnings("ignore")

if len(sys.argv) != 4:
    raise ValueError("Usage: python "+sys.argv[0]+" <working directory> <path to mafft> <path to blast>")

workingdir = sys.argv[1] # working directory path
mafft_path = sys.argv[2]
blast_path = sys.argv[3]
subfamilies = ['L1PA2','L1PA3','L1PA4']

for te in subfamilies:
    # read in the csv with the coordinate and sequence data for each TE
    df = pd.read_csv(workingdir+'/sequences/'+te+'.csv', index_col=0, sep='\t')
    columns = df.columns

    with open(workingdir+'/orf/'+te+'.fasta','w') as blast_fasta:
        for index,row in df.iterrows():
            with open(workingdir+'/alignments/intemp.fa','w') as temp_fasta:
                for i in range(len(row)):
                    if i % 2 != 0:
                        w_str = '>' + columns[i-1] + ' ' + row[i-1] + '\n' + row[i] + '\n'
                        num = row[i-1].split('_')[-2]
                        temp_fasta.write(w_str)
                        blast_fasta.write(w_str)
            # here add details for running mafft
            os.system(mafft_path+' --globalpair --quiet --maxiterate 1000 '+workingdir+'/alignments/intemp.fa > '+workingdir+'/alignments/align_'+te+'_'+str(num)+'.fasta')
    os.system(blast_path+'/makeblastdb -in '+workingdir+'/orf/'+te+'.fasta -dbtype nucl -out '+workingdir+'/orf/'+te)
    os.system(blast_path+'/tblastn -query '+workingdir+'/orf/ORFS.fasta -db '+workingdir+'/orf/'+te+' -out '+workingdir+'/orf/blast_results_'+te+'.txt -outfmt 0')

