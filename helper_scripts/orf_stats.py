import pandas as pd

def read_orf_file(path):
    te_list = ['L1PA2','L1PA3','L1PA4']
    data_storage = dict()
    orfs = ['orf1_ref','orf2_ref']

    # initialize the dictionary to store the data
    for te in te_list:
        data_storage[te] = dict()
        for orf in orfs:
            data_storage[te][orf] = pd.DataFrame(columns=['TE label','orf start','orf end',
                                                            'TE start','TE end','identity','translated TE','orf seq'])
        
    for te in te_list:
        orf = ''
        label = ''
        frame = ''
        spec = ''
        orf_start, orf_end = float('inf'), 0
        TE_start, TE_end = float('inf'), 0
        per_ind = 0
        translated_TE = ''
        orf_str = ''
        visited = set()
        pass_block = 0
        start_comp = lambda x, y: x > y
        end_comp = lambda x, y: x < y
        first = 1
        with open(path+'/orf/results_'+te+'.fasta','r') as align_file:
            for line in align_file:
                if line[:6] == 'Query=':
                    # get the orf label
                    orf = line.split()[-1]
                    visited = set()

                elif line[0] == '>':
                    # found the te label
                    if not first:
                        new_row = {'TE label':label,
                                   'species':spec,
                                'orf start':orf_start,
                                'orf end':orf_end,
                                'TE start':TE_start,
                                'TE end':TE_end,
                                'identity':per_ind,
                                'translated TE':translated_TE,
                                'orf seq': orf_str}
                        new_df = pd.DataFrame([new_row])
                        data_storage[te][orf] = pd.concat([data_storage[te][orf], new_df], ignore_index=True)
                    else:
                        first = 0
                    spec, label = line[1:].strip().split()
                    label = label.split(':')[1].split('_',1)[-1][:-2]
                    pass_block = 0
                    translated_TE = ''
                    orf_str = ''
                    TE_start,TE_end = float('inf'), 0
                    orf_start,orf_end = float('inf'), 0
                    

                elif line[:5] == ' Iden':
                    if label in visited:
                        # already added one of these so don't add more
                        pass_block = 1
                    else:
                        # save the identity value
                        id_str = line.strip().split(',')[0]
                        translated_TE = ''
                        num, denom = id_str.split()[2].split('/')
                        per_ind = int(num) / int(denom)
                        visited.add(label)
                
                elif not pass_block:
                    if line[:5] == ' Fram':
                        frame = line[8:].strip()
                        if frame[0] == '+':
                            start_comp = lambda x, y: x > y
                            end_comp = lambda x, y: x < y
                        else:
                            TE_start,TE_end = 0, float('inf')
                            orf_start,orf_end = float('inf'), 0
                            start_comp = lambda x, y: x < y
                            end_comp = lambda x, y: x > y

                    elif line[:5] == 'Sbjct':
                        q_vals = line.split()
                        if start_comp(TE_start,int(q_vals[1])):
                            TE_start = int(q_vals[1])
                        
                        if end_comp(TE_end,int(q_vals[3])):
                            TE_end = int(q_vals[3])
                        translated_TE += q_vals[2]

                    elif line[:7] == 'Query  ':
                        q_vals = line.split()
                        if orf_start > int(q_vals[1]):
                            orf_start = int(q_vals[1])
                        
                        if orf_end < int(q_vals[3]):
                            orf_end = int(q_vals[3])
                        orf_str += q_vals[2]

                elif line[:11] == '  Database:':
                    new_row = {'TE label':label,
                               'species':spec,
                            'orf start':orf_start,
                            'orf end':orf_end,
                            'TE start':TE_start,
                            'TE end':TE_end,
                            'identity':per_ind,
                            'translated TE':translated_TE,
                            'orf seq': orf_str}
                    new_df = pd.DataFrame([new_row])
                    data_storage[te][orf] = pd.concat([data_storage[te][orf], new_df], ignore_index=True)
    
    return data_storage

def orf_stats(orf,row):
    orf_found = 0
    star_excluded = 0
    hyph_excluded= 0
    len_excluded=0
    per_ind_excluded= 0
    no_alignment = 0
    print(row)
    if len(row['translated TE']) == 0:
        no_alignment = 1

    elif orf == 'orf1_ref':
        if (row['orf start'] != 1 or row['orf end'] != 338):
            len_excluded = 1
        elif '*' in row['translated TE']:
            # if the stop codon is in the translated DNA sequence
            star_excluded = 1
        elif '-' in row['orf seq'] or '-' in row['translated TE']:
            hyph_excluded = 1
        else:
            spec = row['TE label'].split('_')[0]
            if spec == 'Human' and row['identity'] < 0.95:
                # percent identity should be >95%
                per_ind_excluded = 1
            elif spec != 'Human' and row['identity'] < 0.85:
                # percent identity should be > 85%
                per_ind_excluded = 1
            else:
                orf_found = 1
        
    elif orf == 'orf2_ref':
        if (row['orf start'] != 1 or row['orf end'] < 1163):
            len_excluded = 1 
        elif '*' in row['translated TE'][:1163]:
            # if the stop codon is in the translated DNA sequence
            star_excluded = 1
        elif '-' in row['orf seq'][:1163] or '-' in row['translated TE'][:1163]:
            hyph_excluded = 1
        else:
            spec = row['TE label'].split('_')[0]
            if spec == 'Human' and row['identity'] < 0.95:
                # percent identity should be >95%
                per_ind_excluded = 1
            elif spec != 'Human' and row['identity'] < 0.85:
                # percent identity should be > 85%
                per_ind_excluded = 1
            else:
                orf_found = 1
    
    return [orf_found,star_excluded,hyph_excluded,len_excluded,per_ind_excluded,no_alignment]


# if __name__=='__main__':
#     data = read_orf_file('/Users/rachelparsons/Downloads/T2T_primates/L1PA_ILS_Analysis')
#     print(data['L1PA2']['orf1_ref'].shape)
#     print(data['L1PA2']['orf1_ref'].head())

#     for index, row in data['L1PA2']['orf1_ref'].iterrows():
#         res = orf_stats('orf1_ref',row)
#         print(res)
