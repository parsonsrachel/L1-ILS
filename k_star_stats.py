import dendropy
from helper_scripts.patterns import id_with_pat
from itertools import combinations
import pandas as pd
import re
import sys
from collections import Counter

def convert_to_quartet(tree_str):
    new_str = tree_str.replace('[&U] ','')
    if new_str.count('(') > 1:
        # not a star tree
        taxa = dendropy.TaxonNamespace()
        tree = dendropy.Tree.get(data=new_str, schema='newick', taxon_namespace=taxa)
        tree.encode_bipartitions()
        ret_str1, ret_str2 = '', ''

        for edge in tree.internal_edges():
            bipart = edge.bipartition
            group2 = [taxon.label for taxon in taxa if bipart.split_bitmask & taxa.taxon_bitmask(taxon)]
            group1 = [taxon.label for taxon in taxa if taxon.label not in group2]
            g1_str = ','.join(sorted(group1))
            g2_str = ','.join(sorted(group2))
            ret_str1 = f"{g1_str}|{g2_str}"
            ret_str2 = f"{g2_str}|{g1_str}"
        return ret_str1, ret_str2
    
    else:
        taxa = dendropy.TaxonNamespace()
        tree = dendropy.Tree.get(data=new_str, schema='newick', taxon_namespace=taxa)
        leaves = [leaf.taxon.label for leaf in list(tree.leaf_node_iter())]
        return ','.join(sorted(leaves)), None
    
def load_data(filename):
    df = pd.read_csv(filename)
    filtered = df[(df['lifted_all'] == 1) & (df['segdup?'] == 0)]
    orths = filtered[filtered['orthologous'].str.contains('0')]
    orth_filter = pd.concat([filtered, orths, orths]).drop_duplicates(keep=False)
    no_pa = orth_filter[~orth_filter['presence/absence'].str.contains(r'\?', regex=True)]

    drop_list = []
    for index, row in no_pa.iterrows():
        drop = 0
        pa = list(map(int, row['presence/absence'].split(',')))
        if 'None' not in row['TSD_pattern']:
            tsd = list(map(int, row['TSD_pattern'].split(','))) 
        else:
            tsd = [0] * 6
            drop = 1
        for i in range(len(pa)):
            if pa[i] == 1 and tsd[i] != 1:
                add = 0
                drop = 1
        if drop:
            drop_list.append(index)

    no_pa = no_pa.drop(drop_list)
    return no_pa

def threshold_tree(treestr,threshold):
    tree = dendropy.Tree.get(data=re.sub(r"/([0-9.]+)", r"\1",treestr), schema="newick")
    tree.encode_bipartitions()

    for node in tree.postorder_node_iter():
        node.edge.length = None
        if node.is_leaf() or node.label is None:
            continue
        try:
            support = float(node.label)
            if support < threshold:
                node.edge.collapse()
        except ValueError:
            pass  # skip nodes with non-numeric labels

    return tree

def geneTrees(df, pat):
    
    taxa = ['B.Orangutan', 'S.Orangutan', 'Gorilla', 'Human', 'Bonobo', 'Chimp']
    print(f'pattern: {pat}')
    # find all quartet topologies
    topologies = dict()

    # Generate all 4-taxon subsets
    quartet_sets = list(combinations(taxa, 4))
    # get all topologies on 4 species
    for quartet in quartet_sets:
        visited = []
        topologies[','.join(sorted(quartet))] = 0 # add the star tree
        for q in list(combinations(quartet,2)):
            if set(q) not in visited:
                temp_q = set(quartet)
                visited.append(set(q))
                # find the complementary pair of species
                for i in q:
                    temp_q.remove(i)
                visited.append(temp_q)
                group1 = ','.join(sorted(q))
                group2 = ','.join(sorted(temp_q))

                # Combine into a  string
                quartet_str = f"{group1}|{group2}"
                topologies[quartet_str] = 0

    # df_align_store[topologies] = np.nan # add columns to df
    new_df = df[df['presence/absence'] == ','.join(str(x) for x in pat)]
    num_trees = 0
    subset = ['B.Orangutan', 'Chimp', 'Gorilla', 'Human']
    subset_labels = dict()
    # print(new_df.head())
    for _, row in new_df.iterrows():
        # print(row['gene_tree'])
        num_trees += 1
        tree = threshold_tree(row['gene_tree'],0.9)
        # tree = dendropy.Tree.get(data=re.sub(r"/([0-9.]+)", r"\1",row['gene_tree']), schema="newick")
        # threshold = 0.9
        # tree.encode_bipartitions()

        # for node in tree.postorder_node_iter():
        #     node.edge.length = None
        #     if node.is_leaf() or node.label is None:
        #         continue
        #     try:
        #         support = float(node.label)
        #         if support < threshold:
        #             node.edge.collapse()
        #     except ValueError:
        #         pass  # skip nodes with non-numeric labels

        for quartet in quartet_sets:
            sub_tree = tree.extract_tree_with_taxa_labels(quartet)
            sub_tree.encode_bipartitions()
            combo1, combo2 = convert_to_quartet(sub_tree.as_string(schema="newick"))
            # if len()
            # side1, side2 = combo1.split('|')
            # l=side1.split(',')+side2.split(',')
            if '|' not in combo1:
                temp_combo = combo1
                threshold = 85
                while '|' not in temp_combo:
                    temp_tree = threshold_tree(row['gene_tree'],threshold/100)
                    temp_sub_tree = temp_tree.extract_tree_with_taxa_labels(quartet)
                    temp_sub_tree.encode_bipartitions()
                    temp_combo, _ = convert_to_quartet(temp_sub_tree.as_string(schema="newick"))
                    if '|' not in temp_combo:
                        threshold -= 5
                # print(threshold)
                s = combo1.split(',')
                s.sort()
                if s == subset:
                    # print('in if')
                    # if f'{combo1}:{threshold}' in subset_labels:
                    #     subset_labels[f'{combo1}:{threshold}'].add(row['TE_label'])
                    # elif f'{combo2}:{threshold}' in subset_labels:
                    #     subset_labels[f'{combo2}:{threshold}'].add(row['TE_label'])
                    # elif combo1 in topologies:
                    #     subset_labels[f'{combo1}:{threshold}'] = set()
                    #     subset_labels[f'{combo1}:{threshold}'].add(row['TE_label'])
                    # elif combo2 in topologies:
                    #     subset_labels[f'{combo2}:{threshold}'] = set()
                    #     subset_labels[f'{combo2}:{threshold}'].add(row['TE_label'])
                    if combo1 in subset_labels:
                        subset_labels[combo1].add(row['TE_label'])
                    elif combo2 in subset_labels:
                        subset_labels[combo2].add(row['TE_label'])
                    elif combo1 in topologies:
                        subset_labels[combo1] = set()
                        subset_labels[combo1].add(row['TE_label'])
                    elif combo2 in topologies:
                        subset_labels[combo2] = set()
                        subset_labels[combo2].add(row['TE_label'])
            elif '|' in combo1:
                side1, side2 = combo1.split('|')
                l=side1.split(',')+side2.split(',')
                l.sort()
                if l == subset:
                    # print('in second if')
                    if combo1 in subset_labels:
                        subset_labels[combo1].add(row['TE_label'])
                    elif combo2 in subset_labels:
                        subset_labels[combo2].add(row['TE_label'])
                    elif combo1 in topologies:
                        subset_labels[combo1] = set()
                        subset_labels[combo1].add(row['TE_label'])
                    elif combo2 in topologies:
                        subset_labels[combo2] = set()
                        subset_labels[combo2].add(row['TE_label'])


            if combo1 in topologies:
                topologies[combo1] +=1

            elif combo2 in topologies: 
                topologies[combo2] += 1
        
    # topologies contains the quartets and their counts
    # print(unresolved)
    # assert(len(unresolved)==1547)
    top_df = pd.DataFrame(list(topologies.items()), columns=["Quartet", "Count"])
    return subset_labels, top_df

            

if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError("Usage: python "+sys.argv[0]+" <working directory>")

    path = sys.argv[1]
    df2 = load_data(path+'/call_set/L1PA2.csv')
    df3 = load_data(path+'/call_set/L1PA3.csv')
    df4 = load_data(path+'/call_set/L1PA4.csv')
    new_df = pd.DataFrame(columns = ['subset','label','te_region_length','all_non_gap','all_same','singleton','quartet','triplet','all_unique',
                                     'percent all the same','percent singleton','percent quartet','percent triplet','percent all unique', 'percent uninformative', 'percent informative'])

    all_df = pd.concat([df3, df4], ignore_index=True)
    four_df = pd.concat([df2, df3, df4], ignore_index=True)
    all_labels, _ = geneTrees(all_df,(1,1,1,1,1,1))
    # four_unresolved, _ = geneTrees(four_df,(0,0,1,1,1,1))
    # print(all_unresolved)
    for quart in all_labels:
        for label in all_labels[quart]:
            # print(label)
            if label[4] == '2':
                df = df2
            elif label[4] == '3':
                df = df3
            elif label[4] == '4':
                df = df4
            r_len = df.loc[df['TE_label'] == label, 'region_length'].values[0]

            # msa = {k: '' for k in quart.split(',')}
            msa = dict()
            with open(f'{path}/alignments/align_{label}.fasta','r') as alignment:
                save = 0
                spec = ''
                for line in alignment:
                    if line[0] == '>':
                        spec = line.split()[0][1:]
                        msa[spec] = ''
                    else:
                        msa[spec] += line.strip()
            # print(msa)
            start_index, end_index = 0, len(msa['Human'])-1
            count = 0
            for c in range(len(msa['Human'])):
                if msa['Human'][c] != '-':
                    count+=1
                if count == 500:
                    start_index = c

            count= 0
            for c in range(len(msa['Human'])-1,0,-1):
                if msa['Human'][c] != '-':
                    count+=1
                if count == 500:
                    end_index = c
            #if label == 'L1PA4_2456':
                #print(f'length of human sequence: {len(msa["Human"])}')
                #print(f'start index: {start_index}')
                #print(f'end index: {end_index}')
            all_non_gap = 0
            all_same = 0
            singleton = 0
            quartet = 0
            three_states = 0
            all_unique = 0
            for i in range(start_index+1,end_index,1):
                if '|' in quart:
                    part1, part2 = quart.split('|')
                    species = part1.split(',') + part2.split(',')
                else:
                    species = quart.split(':')[0].split(',')
                chars = [msa[species[0]][i].upper(),msa[species[1]][i].upper(),msa[species[2]][i].upper(),msa[species[3]][i].upper()]
                if chars[0] != '-' and chars[1] != '-' and chars[2] != '-' and chars[3] != '-':
                    all_non_gap += 1
                    counts = sorted(Counter(chars).values())
                    if counts == [2,2]:
                        quartet += 1
                    elif counts == [1,3]:
                        singleton += 1
                    elif counts == [1, 1, 2]:
                        three_states += 1
                    elif counts == [4]:
                        all_same += 1
                    else:
                        all_unique += 1

            new_row = pd.DataFrame([{'subset': quart, "label": label, "te_region_length": r_len, 'all_non_gap': all_non_gap, 'all_same': all_same,
                                      'singleton': singleton, 'quartet': quartet, 'triplet': three_states, 'all_unique': all_unique, 
                                      'percent all the same': all_same/all_non_gap,'percent singleton': singleton/all_non_gap,
                                      'percent quartet': quartet/all_non_gap,'percent triplet':three_states/all_non_gap,'percent all unique':all_unique/all_non_gap,
                                      'percent uninformative': (all_same+singleton+all_unique)/all_non_gap, 'percent informative': (quartet+three_states)/all_non_gap}])
                                      
            new_df = pd.concat([new_df, new_row], ignore_index=True)

    new_df.to_csv(f'{path}/figures/star_stats.csv')


