import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from helper_scripts.filterTEs import filter
import argparse

parser = argparse.ArgumentParser(description="Isolate the TE annotations for a specific TE subfamily and write the data to a bed file.")

# Add arguments
parser.add_argument("-i","--input", type=str,
                    help="Path to the input bigBed file",required=True)
parser.add_argument("-o1", "--output1", type=str,
                    help="Path to the directory of liftover TE coordinates", required=True)
parser.add_argument("-o2", "--output2", type=str,
                    help="Path to the directory of filtered RM data callset", required=True)
parser.add_argument("-r", "--repeat", type=str, nargs="+",
                     help="Subfamily label of the annotations to isolate")

# Parse arguments
args = parser.parse_args()

# turn big bed file into repeatmasker out format
rm_format_file = args.input[:args.input.rfind('/', 0, len(args.input))+1]+'rm_annotations.out'
with open(args.input,'r') as bed, open(rm_format_file,'w') as rm:
    for line in bed:
        lst = line.split()[13:]
        for l in lst:
            if ',' in l:
                rm.write('\n'+l[1:]+'\t')
            else:
                rm.write(l+'\t')
        rm.write('\n')

# read in the repeatmasker results for chm13 and filter
df = pd.read_csv(rm_format_file,
                sep=r'\s+',  
                names=[
                    "SW_score", "perc_div", "perc_del", "perc_ins",
                    "query_sequence", "query_begin", "query_end", "query_left",
                    "strand", "repeat_name", "repeat_class_family",
                    "repeat_begin", "repeat_end", "repeat_left", "ID","star"
                ])

# remove the sex chromosomes
df = df[df["query_sequence"]!='chrX']
df = df[df["query_sequence"]!='chrY']


te_df = df[(df["repeat_name"]==args.repeat[0]) | (df["repeat_name"]==args.repeat[1]) | (df["repeat_name"]==args.repeat[2])]
TE, total, filt1, filt2 = filter(te_df) # filter based on sequence length and distance to closest TE
print(f'in total there are {total[0]}, {total[1]} ,{total[2]}  L1PA2-4 elements')
print(f'{total[0] - filt1[0]}, {total[1] - filt1[1]}, {total[2] - filt1[2]} L1 elements removed by first filter (len >= 200)')
print(f'{filt1[0] - filt2[0]}, {filt1[1] - filt2[1]}, {filt1[2] - filt2[2]} L1 elements removed after second filter (no element within 500bp)')

for te in args.repeat:
    single_te = TE[TE["repeat_name"]==te].reset_index(drop=True)
    single_te['TE_label'] = [f'{te}_{i}' for i in single_te.index]
    print(f'{single_te.shape[0]} {te} elements remaining')
    single_te.to_csv(args.output2+ '/' + te + '.csv', index = False)

    bed_df = single_te[['query_sequence','query_begin','query_end']]
    bed_df['query_begin'] = bed_df['query_begin'] - 500
    bed_df['query_end'] = bed_df['query_end'] + 500
    bed_df['query_sequence'] = 'chm13#1#' + bed_df['query_sequence']
    bed_df['name'] = [f'{te}_{i}' for i in single_te.index]
    bed_df[['ID']] = single_te[['ID']]
    bed_df['strand'] = '+'

    bed_df.to_csv(args.output1 + '/' + te + '.bed', index=False, header=False,  sep='\t') 

