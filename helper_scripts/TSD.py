
def get_n_mer(seq_str, flank, n):
    '''
    Input: 
    species = label of the speces we want to find the TSD for
    msa = the output of the getMSA function {species: (coordinates, sequence)} for the range we want to find all n-mers for
    n = size of the kmer
    
    Output:
    list of all nmers
    '''
    seq = seq_str.replace('-','')
    if flank == 1:
        seq = seq[-65:].upper()
    elif flank == 2:
        seq = seq[:65].upper()
    nmers = set()
    for i in range(0,(len(seq)-n)+1):
        if i == len(seq)-n:
            nmers.add(seq[i:])
        else:
            nmers.add(seq[i:i+n])
    return nmers

def per_AT(seq):
    '''
    Input: sequence
    Output: percent of A/T in the seq
    '''

    a_count, t_count = 0,0
    for i in seq:
        if i.upper() == 'A':
            a_count += 1
        if i.upper() == 'T':
            t_count += 1
    if a_count + t_count == len(seq):
        return 1
    elif a_count > t_count:
        return a_count / len(seq)
    else:
        return t_count / len(seq)

def tsd_human(tsd1,tsd2):
    for n in range(30,4,-1):
        lead_mers = get_n_mer(tsd1, 1, n)
        tail_mers = get_n_mer(tsd2, 2, n)
        if lead_mers.intersection(tail_mers):
            mer_set = lead_mers & tail_mers
            while len(mer_set) != 0:
                mer_str = mer_set.pop()
                if per_AT(mer_str) < 0.8:
                    lead_index_offset = tsd1.upper().find(mer_str)
                    tail_index_offset = tsd2.upper().find(mer_str)

                    return mer_str, lead_index_offset, tail_index_offset
    return None, None, None

def tsd_non_human(msa,tsd1_start,tsd2_start,seq): # start and stop ar the indices where the tsd is in the human (with no gaps)
    tsd_len = len(seq)
    nuc_count = 0
    msa_t1 = [0,0]
    msa_t2 = [0,0]
    

    for i in range(len(msa['Human'])):
        # find just the ailgnment in the tsd region
        if msa['Human'][i] != '-':
            nuc_count+=1

        if nuc_count == tsd1_start + 1:
            msa_t1[0] = i

        elif nuc_count == tsd1_start + 1 + tsd_len and msa_t1[1] == 0:
            msa_t1[1] = i

        elif nuc_count == tsd2_start + 1:
            msa_t2[0] = i

        elif nuc_count == tsd2_start + 1 + tsd_len:
            msa_t2[1] = i
            break
    
    ret_dict = dict()
    h = msa['Human']
    h1 = h[msa_t1[0]:msa_t1[1]]
    h2 = h[msa_t2[0]:msa_t2[1]]
    for spec in msa:
        if spec == 'Human':
            ret_dict[spec] = (h1,h2,1)
        else:
            non_gap_t1, non_gap_t2 = 0,0
            match_t1, match_t2 = 0,0
            nuc_hum_t1, nuc_hum_t2 = 0,0
            # print(spec)
            s1 = msa[spec][msa_t1[0]:msa_t1[1]]
            s2 = msa[spec][msa_t2[0]:msa_t2[1]]
            # print(f'{s1} to {h1}')
            # print(f'{s2} to {h2}')

            for ind in range(len(s1)):
                if s1[ind].upper() == h1[ind].upper():
                    match_t1 += 1
                if s1[ind] != '-' and h1[ind] != '-':
                    non_gap_t1 += 1
                if h1[ind] != '-':
                    nuc_hum_t1 += 1
            
            for ind in range(len(s2)):
                if s2[ind].upper() == h2[ind].upper():
                    match_t2 += 1
                if s2[ind] != '-' and h2[ind] != '-':
                    non_gap_t2 += 1
                if h2[ind] != '-':
                    nuc_hum_t2 += 1
            # print(f'matches 1: {match_t1} 2: {match_t2}')
            # print(f'non_gap 1: {non_gap_t1} 2: {non_gap_t2}')
            # print(f'len of region 1: {msa_t1[1] - msa_t1[0]} 2: {msa_t2[1] - msa_t2[0]}')
            if non_gap_t1 / nuc_hum_t1 > 0.8 and non_gap_t2 / nuc_hum_t2 > 0.8 and match_t1 / non_gap_t1 > 0.75 and match_t2 / non_gap_t2 > 0.75:
                ret_dict[spec] = (s1,s2,1)
            else: 
                ret_dict[spec] = (s1,s2,0)

    return ret_dict

# if __name__=='__main__':

#     msa = dict()
#     with open('/Users/rachelparsons/Downloads/T2T_primates/L1PA_ILS_Analysis/alignments/align_L1PA2_182.fasta','r') as fa:
#         spec = ''
#         start, stop = 0, 0
#         for line in fa:
#             if line[0] == '>':
#                 spec = line.split()[0][1:]
#                 if spec == 'Human':
#                     # get the start and stops
#                     start, stop = line.split()[1].split(':')[1].split('_')[0].split('-')
#             else:
#                 if spec in msa:
#                     msa[spec] += line.strip()
#                 else:
#                     msa[spec] = line.strip()
    
#     # print(msa)

#     tsd1_region = msa['Human'].replace('-','')[450:515]
#     tsd2_region = msa['Human'].replace('-','')[-515:-450]
#     # print(tsd1_region)
#     # print(tsd2_region)

#     seq, lead_off, tail_off = tsd_human(tsd1_region,tsd2_region)
#     print(seq,lead_off,tail_off)

#     res = tsd_non_human(msa,lead_off+450,(-515+tail_off)+len(msa['Human'].replace('-','')),seq)
#     for s in res:
#         print(s)
#         print(res[s])

