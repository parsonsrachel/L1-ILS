import sys
import pandas as pd
from helper_scripts.filterTEs import *

def no_mapping(row):
    cols = ['S.Orangutan','B.Orangutan','Gorilla_merged','Bonobo_merged','Chimp']
    for col in cols:
        try:
            if float(row[col]) < 0:
                return True
        except (ValueError, TypeError):
            continue
    return False

if __name__=='__main__':
    if len(sys.argv) != 3:
        raise ValueError("Usage: python "+sys.argv[0]+" <working directory> <seg dup directory>")

    path = sys.argv[1] # working directory path
    seg_dup_path = sys.argv[2]

    te_list = ['L1PA2','L1PA3','L1PA4']
    # load all the segmental duplication data
    cols = ['Chimp','S.Orangutan','B.Orangutan','Bonobo_merged','Gorilla_merged']
    long_cols = ['query_sequence','query_begin','query_end','Chimp','S.Orangutan','B.Orangutan','Bonobo_merged','Gorilla_merged']
    labels = ['Human','Chimp','S.Orangutan','B.Orangutan','Bonobo','Gorilla'] 

    segs = {seg_dup_path+'/AG05252_PPY.SDs.bed':'B.Orangutan',
            seg_dup_path+'/AG06213_PAB.SDs.bed':'S.Orangutan',
            seg_dup_path+'/AG18354_PTR.SDs.bed':'Chimp',
            seg_dup_path+'/HG002.SDs.bed':'Human',
            seg_dup_path+'/Jim_GGO.SDs.bed':'Gorilla',
            seg_dup_path+'/PR00251_PPA.SDs.bed':'Bonobo'}

    # read in segmental duplications
    small_cols = ['chr','begin','end']
    sd_dict = {'B.Orangutan':pd.DataFrame(columns=small_cols),
                'S.Orangutan':pd.DataFrame(columns=small_cols),
                'Chimp':pd.DataFrame(columns=small_cols),
                'Human':pd.DataFrame(columns=small_cols),
                'Gorilla':pd.DataFrame(columns=small_cols),
                'Bonobo':pd.DataFrame(columns=small_cols)}

    for s_file in segs:
        with open(s_file,'r') as sd_file:
            for line in sd_file:
                sp = line.split()
                sd_dict[segs[s_file]].loc[len(sd_dict[segs[s_file]])] = {'chr': sp[0], 'begin': int(sp[1]), 'end': int(sp[2])} 

    # loop over the unmapped first and mark in table as 
    unmap_dict = {'#Deleted in new':-1,
                '#Duplicated in new':-2,
                '#Partially deleted in new':-3,
                '#Split in new':-4}


    for te in te_list:
        df = pd.read_csv(path+'/call_set/'+te+'.csv')

        # mark all the times liftover failed to find coordinates
        unmapped = {path+'/lifted_beds/'+te+'/chm13-to-mPonAbe1#1-unmapped.bed':'S.Orangutan',
                    path+'/lifted_beds/'+te+'/chm13-to-mPonPyg2#1-unmapped.bed':'B.Orangutan',
                    path+'/lifted_beds/'+te+'/chm13-to-mGorGor1#M-unmapped.bed':'Gorilla M',
                    path+'/lifted_beds/'+te+'/chm13-to-mGorGor1#P-unmapped.bed':'Gorilla P',
                    path+'/lifted_beds/'+te+'/chm13-to-mPanPan1#M-unmapped.bed':'Bonobo M',
                    path+'/lifted_beds/'+te+'/chm13-to-mPanPan1#P-unmapped.bed':'Bonobo P',
                    path+'/lifted_beds/'+te+'/chm13-to-mPanTro3#1-unmapped.bed':'Chimp'}
        for filename in unmapped:
            # print(filename)
            df[unmapped[filename]] = 1
            with open(filename,'r') as unmap_file:
                prev_line = ''
                for line in unmap_file:
                    if line[0] == '#':
                        prev_line = line.strip()
                    else:
                        parts = line.strip().split()
                        id, label = int(parts[4]), parts[3]
                        df.loc[(df['TE_label'] == label), unmapped[filename]] = unmap_dict[prev_line]
 
        mapped = {path+'/lifted_beds/'+te+'/chm13-to-mPonAbe1#1.bed':'S.Orangutan',
                    path+'/lifted_beds/'+te+'/chm13-to-mPonPyg2#1.bed':'B.Orangutan',
                    path+'/lifted_beds/'+te+'/chm13-to-mGorGor1#M.bed':'Gorilla M',
                    path+'/lifted_beds/'+te+'/chm13-to-mGorGor1#P.bed':'Gorilla P',
                    path+'/lifted_beds/'+te+'/chm13-to-mPanPan1#M.bed':'Bonobo M',
                    path+'/lifted_beds/'+te+'/chm13-to-mPanPan1#P.bed':'Bonobo P',
                    path+'/lifted_beds/'+te+'/chm13-to-mPanTro3#1.bed':'Chimp'}
    
        for s in mapped:
            with open(s,'r') as bed:
                for line in bed:
                    sp = line.split()
                    new_str = sp[0] + ':' + sp[1] + '-' + sp[2] + '_' + sp[3] + '_' + sp[4] + '_' + sp[5]
                    df.loc[(df['TE_label'] == sp[3]), mapped[s]] = new_str

 
        # get the haplotypes for the gorilla and bonobo
        bonobo =['chr1_pat_hsa1','chr12_mat_hsa2a','chr13_mat_hsa2b','chr2_mat_hsa3',
                'chr3_pat_hsa4','chr4_mat_hsa5','chr5_pat_hsa6','chr6_mat_hsa7',
                'chr7_pat_hsa8','chr11_mat_hsa9','chr8_mat_hsa10','chr9_pat_hsa11',
                'chr10_mat_hsa12','chr14_pat_hsa13','chr15_mat_hsa14','chr16_pat_hsa15',
                'chr18_pat_hsa16','chr19_pat_hsa17','chr17_mat_hsa18','chr20_pat_hsa19',
                'chr21_mat_hsa20','chr22_pat_hsa21','chr23_pat_hsa22']
        gorilla=['chr1_pat_hsa1','chr11_mat_hsa2b','chr12_pat_hsa2a','chr2_pat_hsa3',
                'chr3_pat_hsa4','chr19_pat_hsa5x17','chr4_pat_hsa17x5','chr5_mat_hsa6',
                'chr6_mat_hsa7','chr7_pat_hsa8','chr13_pat_hsa9','chr8_pat_hsa10',
                'chr9_pat_hsa11','chr10_mat_hsa12','chr14_pat_hsa13','chr15_pat_hsa14',
                'chr16_pat_hsa15','chr18_pat_hsa16','chr19_pat_hsa5x17','chr4_pat_hsa17x5',
                'chr17_mat_hsa18','chr20_mat_hsa19','chr21_pat_hsa20','chr22_mat_hsa21','chr23_mat_hsa22']

        index_mat = df['Bonobo M'].apply(
                        lambda x: x.split('#')[-1].split(':')[0] if isinstance(x, str) else x
                    ).isin(bonobo)

        index_pat = df['Bonobo P'].apply(
                        lambda x: x.split('#')[-1].split(':')[0] if isinstance(x, str) else x
                    ).isin(bonobo)

        df['Bonobo_merged'] = df.apply(
            lambda row: row['Bonobo M'] if index_mat[row.name]
            else row['Bonobo P'] if index_pat[row.name]
            else row['Bonobo M'] if isinstance(row['Bonobo M'], int)
            else row['Bonobo P'],
            axis=1
        )
        df = df.drop('Bonobo M', axis=1)
        df = df.drop('Bonobo P', axis=1)

        index_mat = df['Gorilla M'].apply(
                        lambda x: x.split('#')[-1].split(':')[0] if isinstance(x, str) else x
                    ).isin(gorilla)

        index_pat = df['Gorilla P'].apply(
                        lambda x: x.split('#')[-1].split(':')[0] if isinstance(x, str) else x
                    ).isin(gorilla)

        df['Gorilla_merged'] = df.apply(
            lambda row: row['Gorilla M'] if index_mat[row.name]
            else row['Gorilla P'] if index_pat[row.name]
            else row['Gorilla M'] if isinstance(row['Gorilla M'], int)
            else row['Gorilla P'],
            axis=1
        )

        df = df.drop('Gorilla M', axis=1)
        df = df.drop('Gorilla P', axis=1)
        df['lifted_all'] = df.apply(lambda row: 1 if not no_mapping(row) else 0, axis=1)
        count_lifted_all = df['lifted_all'].sum()
        print(f'{te} has {count_lifted_all} elements that mapped, {df.shape[0] - count_lifted_all} that failed to map, {df.shape[0]} total elements')
 

        # label all the elements in segmental duplications
        df['segdup?'] = 0
        indeces = df[cols].apply(
                lambda row: all(isinstance(val, str) for val in row),
                axis=1
            )
        filtered = df[long_cols][indeces]
        # print(filtered)
        for s in labels:
            col_str = ''
            if s == 'Human':
                rows = filtered[['query_sequence', 'query_begin', 'query_end']]
            elif s == 'Gorilla' or s == 'Bonobo':
                col_str = s+'_merged'
                rows = filtered[[col_str]]
            else:
                col_str = s
                rows = filtered[[col_str]]
            
            for ind,r in rows.iterrows():
                if s == 'Human':
                    # make sure the flanking region is included for human
                    chr = r['query_sequence']
                    start = str(int(r['query_begin']) - 500)
                    stop = str(int(r['query_end']) + 500)
                else:
                    val = r[col_str]
                    chr = val.split('#')[-1].split(':')[0]
                    start, stop = val.split('#')[-1].split(':')[1].split('_')[0].split('-')

                mask = ((sd_dict[s]['chr'] == chr) & 
                        (sd_dict[s]['begin'] <= int(stop)) & 
                        (sd_dict[s]['end'] >= int(start)))
                is_segdup = mask.any()

                if is_segdup:
                    df.at[ind, 'segdup?'] = 1
        count_in_segdup = df['segdup?'].sum()
        print(f'{te} had {count_in_segdup} in segmental duplication region')

        n_cols = ['Chimp','S.Orangutan','B.Orangutan','Bonobo_merged','Gorilla_merged','lifted_all','segdup?']
        indeces = df[n_cols].apply(
                lambda row: all(isinstance(val, str) and row[-2] == 1 and row[-1] == 0 for val in row[:-2]),
                axis=1
            )
        filtered = df[indeces]
        assert((filtered['segdup?'] == 1).any() == False)

        for c in cols:
            new_c = c.split('_merged')[0]
            filtered[[c]].to_csv(path+'/coords_to_pull/'+te+'/'+new_c+'.csv', index=False, header=False, sep=' ')
        human = filtered[['query_sequence','query_begin','query_end','TE_label']]
        # adjust the human coordinates to include the flankers
        human['query_begin'] = human['query_begin'] - 500
        human['query_end'] = human['query_end'] + 500

        human.to_csv(path+'/coords_to_pull/'+te+'/Human.csv', index=False, header=False, sep=' ')

        df.to_csv(path+'/call_set/'+te+'.csv',index=False)

    

