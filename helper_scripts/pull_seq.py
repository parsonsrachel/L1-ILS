def rv_comp(s):
    rev_s = s[::-1]
    rev_comp_s = ''
    for i in rev_s.upper():
        if i == 'A':
            rev_comp_s+= 'T'
        elif i == 'T':
            rev_comp_s+= 'A'
        elif i == 'C':
            rev_comp_s+= 'G'
        else:
            rev_comp_s+= 'C'
    return rev_comp_s

def pull(genome_file):
    read_genome = dict()
    with open(genome_file,'r') as genome:
        chr_num = ''
        for line in genome:
            line = line.strip()
            if line[0] == '>':
                # find which chromosome
                chr_num = line[1:]
                read_genome[chr_num] = []
            else:
                read_genome[chr_num].append(line)
    # Join sequences together
    for chr_num in read_genome:
        read_genome[chr_num] = ''.join(read_genome[chr_num])

    return read_genome
