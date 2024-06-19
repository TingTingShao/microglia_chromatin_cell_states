#!/usr/bin/env python

import os, sys, argparse
import time
import pandas as pd
import pickle

start_time = time.time()

input_filename = None
output_filenames = None

parser=argparse.ArgumentParser(description="export gene cell matrix to tsv", \
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-a", '--input_file', nargs=1, default=None, help="input the anndataset")
parser.add_argument("-o", '--output_files', nargs='+', default=None, help="output file names")
args=parser.parse_args()
input_filename = args.input_file[0]
output_filenames = args.output_files

# load
gene_act = pickle.load(open(input_filename, 'rb'))

# adjust content
df = pd.DataFrame(data=gene_act.mtx, index=gene_act.feature_names, columns=gene_act.cell_names)
sampled_df = df.sample(n=3000, random_state=42)

# write
df.to_csv(output_filenames[0], sep='\t')
sampled_df.to_csv(output_filenames[1], sep='\t')

end_time_read = time.time()
print("Time to process: "+str(end_time_read - start_time)+" seconds")



