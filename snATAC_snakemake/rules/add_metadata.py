#!/usr/bin/env python

import os, sys, argparse
import time
import pandas as pd
import anndata as ad
import pickle 
from scipy.sparse import csr_matrix

start_time = time.time()

input_filenames = None
output_filename = None

parser=argparse.ArgumentParser(description="add metadata", \
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-a", '--input_files', nargs='+', default=None, help="input the filenames")
parser.add_argument("-o", '--output_file', nargs=1, default=None, help="output file name")

args=parser.parse_args()
input_filenames = args.input_files
output_filename = args.output_file[0]

# cell_gene.h5ad
gene_act=ad.read_h5ad(input_filenames[0])

# microglia h5ad file
microglia=ad.read_h5ad(input_filenames[1])

dat=microglia.obs.copy()
new_index = []
for index in dat.index:
    sample, cell = index.split(':')
    new_index.append(f"{cell}___{sample}")
dat.index = new_index
print(gene_act.obs.index.equals(dat.index))

dat=dat.reindex(gene_act.obs.index)
print(gene_act.obs.index.equals(dat.index))
gene_act=ad.AnnData(gene_act.X, obs=dat, var=gene_act.var)
gene_act.write(output_filename, compression="gzip")

end_time_read = time.time()
print("Time to process: "+str(end_time_read - start_time)+" seconds")



