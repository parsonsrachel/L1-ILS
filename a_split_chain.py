import sys

if len(sys.argv) != 2:
    raise ValueError("Usage: python "+sys.argv[0]+" <working directory>")

path = sys.argv[1]
d = dict()

with open(path+'/chains/all_species/chm13#1.p70.aln.chain','r') as chain:
    cur_s = ''
    for line in chain:
        if line[:5] == 'chain':
            # check to see what species
            parts = line.split('\t')[7].split('#')
            temp_s = parts[0]+'#'+parts[1]
            if temp_s in d:
                d[temp_s].append(line)

            else:
                d[temp_s] = [line]

            cur_s = temp_s
        else:
            d[cur_s].append(line)
 
for k in d:
    with open(path+'/chains/by_spec_haplotype/chm13-to-'+k+'.aln.chain','w') as new_file:
        for l in d[k]:
            new_file.write(l)

