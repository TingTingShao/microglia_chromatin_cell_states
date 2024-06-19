#!/usr/bin/env python
import duckdb as db
import pickle
import os, sys, argparse
import pandas as pd
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix
# import snapatac2 as snap 
import scanpy as sc
from scipy import io
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)

from pycisTopic.clust_vis import plot_imputed_features
import time

start_time = time.time()

input_filenames = None
output_dir=None

parser=argparse.ArgumentParser(description="add metadata", \
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-a", '--input_files', nargs="+", default=None, help="input files")
parser.add_argument("-f", '--output_folder', nargs=1, default=None, help="output folder")

args=parser.parse_args()

input_filenames = args.input_files
output_dir=args.output_folder[0]
os.makedirs(output_dir, exist_ok=True)

cell_gene2=ad.read_h5ad(input_filenames[0])
nicol_genes=pd.read_csv(input_filenames[1], index_col=0)
cistopic_obj = pickle.load(open(input_filenames[2], 'rb'))
gene_act = pickle.load(open(input_filenames[3], 'rb'))

all_intersect_variable=set(nicol_genes['gene']).intersection(set(cell_gene2.var_names))
nicol_genes = nicol_genes[nicol_genes['gene'].isin(all_intersect_variable)]

top10_nicol_genes = {}
for x in nicol_genes.index.unique():
    top10_nicol_genes[x]=nicol_genes[nicol_genes.index==x].head(10)['gene'].tolist()
print(top10_nicol_genes)

for celltype in top10_nicol_genes.keys():
    plot_imputed_features(
    cistopic_obj,
    reduction_name='SNAP',
    imputed_data=gene_act,
    features=top10_nicol_genes[celltype], 
    scale=True,
    num_columns=4,
    save=f'{output_dir}/{celltype.replace(" ", "_")}_umap.png'
    )


end_time_read = time.time()
print("Time to process: "+str(end_time_read - start_time)+" seconds")
