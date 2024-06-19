#!/usr/bin/env python

import os, sys, argparse
import time
import pandas as pd
import anndata as ad
import pickle 
from scipy.sparse import csr_matrix
import numpy as np

start_time = time.time()

input_filename = None
output_filename = None

parser=argparse.ArgumentParser(description="export gene cell matrix to h5ad", \
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-a", '--input_file', nargs=1, default=None, help="input the gene_act pickle")
parser.add_argument("-o", '--output_file', nargs=1, default=None, help="output anndata file name")
args=parser.parse_args()
input_filename = args.input_file[0]
output_filename = args.output_file[0]

gene_act = pickle.load(open(input_filename, 'rb'))

feature_cell=csr_matrix(gene_act.mtx, dtype=np.float32)
cell_feature=feature_cell.transpose()

adata = ad.AnnData(cell_feature)
adata.obs_names=gene_act.cell_names
adata.var_names=gene_act.feature_names

adata.write(output_filename, compression="gzip")

end_time_read = time.time()
print("Time to process: "+str(end_time_read - start_time)+" seconds")



