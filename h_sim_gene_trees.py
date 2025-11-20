import os
import pandas as pd
import sys

if len(sys.argv) != 3:
    raise ValueError("Usage: python "+sys.argv[0]+" <working directory> <path to iqtree>")

workingdir = sys.argv[1] # working directory path
iqtree = sys.argv[2]
columns=['B.Orangutan','S.Orangutan','Gorilla','Human','Chimp','Bonobo']

#path = '/Users/rachelparsons/Downloads/T2T primates/L1PA ILS Analysis/alignments'
for l in ['L1PA2','L1PA3','L1PA4']:
    # df = pd.read_csv('../gene_tree_tes/to_tree_'+l+'.csv', header=None)
    # for index, row in df.iterrows():
    df = pd.read_csv(workingdir+'/call_set/'+l+'.csv')
    df['gene_tree'] = None
    for filename in os.listdir(workingdir+'/alignments'):
        msa = dict()
        coords = dict()
        print(filename)
        if filename[6:11] == l:
            label = filename.split('.')[0].split('_',1)[1]
            with open(workingdir+'/alignments/'+filename,'r') as fasta_file:
                spec = ''
                for line in fasta_file:
                    if line[0] == '>':
                        # save the species label
                        sp = line.split()
                        spec, dets = sp[0][1:], sp[1]
                        msa[spec] = ''
                        coords[spec] = dets
                    else:
                        msa[spec] += line.strip()
            
            print(label)
            # print(df.columns)
            # if len(df.loc[df['TE_label'] == label, 'presence/absence']) != 1:
            #     pa = ['0','0','0','0','0','0']
            # else:
            #     pa = df.loc[df['TE_label'] == label, 'presence/absence'].split(',')
            pa = df.loc[df['TE_label'] == label, 'presence/absence'].iloc[0].split(',')
            print(pa)
            # pa = df.loc[df['TE_label'] == label, 'presence/absence'].split(',')
            # isolate te
            num_nucs = 0
            start_ind, end_ind = 0, len(msa['Human'])-1
            for c in range(len(msa['Human'])):
                if num_nucs == 500:
                    start_ind = c
                    break
                
                if msa['Human'][c] != '-':
                    num_nucs += 1

            num_nucs = 0
            for c in range(len(msa['Human'])-1, 0, -1):
                if num_nucs == 500:
                    end_ind = c
                    break
                
                if msa['Human'][c] != '-':
                    num_nucs += 1
            
            assert(len(msa['Human'][:start_ind].replace('-','')) == 500)
            assert(len(msa['Human'][end_ind + 1:].replace('-','')) == 500)

            # write just the TE region to a temp file
            with open(workingdir+'/trees/te/'+filename,'w') as temp_file:
                for index in range(len(columns)):
                    if pa[index] == '1':
                        key = columns[index]
                        if msa[key][start_ind:end_ind].count('-') != len(msa[key][start_ind:end_ind]):
                            temp_file.write('>'+key+' '+coords[key]+'\n'+msa[key][start_ind:end_ind]+'\n')

            # write just the TE + f1 + f2 region to a file
            with open(workingdir+'/trees/te_f1_f2/'+filename,'w') as temp_file:
                for key in msa:
                    temp_file.write('>'+key+' '+coords[key]+'\n'+msa[key]+'\n')

            # write just f1
            with open(workingdir+'/trees/f1/'+filename,'w') as temp_file:
                for key in msa:
                    temp_file.write('>'+key+' '+coords[key]+'\n'+msa[key][:start_ind]+'\n')

            # write f2
            with open(workingdir+'/trees/f2/'+filename,'w') as temp_file:
                for key in msa:
                    temp_file.write('>'+key+' '+coords[key]+'\n'+msa[key][end_ind:]+'\n')
            

            # call IQTree
            os.system(iqtree+' -s '+workingdir+'/trees/te/'+filename+' -m GTR+G --abayes')
            os.system(iqtree+' -s '+workingdir+'/trees/te_f1_f2/'+filename+' -m GTR+G --abayes')
            os.system(iqtree+' -s '+workingdir+'/trees/f1/'+filename+' -m GTR+G --abayes')
            os.system(iqtree+' -s '+workingdir+'/trees/f2/'+filename+' -m GTR+G --abayes')
