#!/usr/bin/env python

"""
Predicted values (y_pred): [1, 1, 1, 0, 1]
True values (y_true):      [1, 1, 0, 0, 1]
Jaccard Index: 0.75 3/4
Jaccard Distance: 1/4 > 0.25
https://www.geeksforgeeks.org/how-to-calculate-jaccard-similarity-in-python/

/data/leuven/351/vsc35107/master_thesis/workflow_new/rules/jaccard.py 
-a /lustre1/project/stg_00079/students/tingting/data/sun/snap2_allfragments/08/030results_cistopic/cell_gene02.h5ad 
-o /lustre1/project/stg_00079/students/tingting/data/sun/snap2_allfragments/08/exploratory/jaccard_snap_pycis.png

./jaccard.py 
-a /lustre1/project/stg_00079/students/tingting/data/sun/snap2_allfragments/08/030results_cistopic/cell_gene04.h5ad 
-o /lustre1/project/stg_00079/students/tingting/data/sun/snap2_allfragments/08/exploratory/jaccard_snap_pycis.png

./jaccard.py -a /lustre1/project/stg_00079/students/tingting/data/sun/snap2_allfragments/08/030results_cistopic/cell_gene04.h5ad -o /lustre1/project/stg_00079/students/tingting/data/sun/snap2_allfragments/08/exploratory/jaccard_act_act.png
"""  

import os, sys, argparse
import time
import pandas as pd
import anndata as ad
import pickle 
from scipy.sparse import csr_matrix
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns

start_time = time.time()

input_filename = None
output_filename = None

parser=argparse.ArgumentParser(description="add metadata", \
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-a", '--input_file', nargs=1, default=None, help="input the gene act h5ad (cell_gene.h5ad)")
parser.add_argument("-o", '--output_file', nargs=1, default=None, help="output file name output image name")

args=parser.parse_args()
input_filename = args.input_file[0]
output_filename = args.output_file[0]


cell_gene=ad.read_h5ad(input_filename)

#TODO hard coding parameters!
df_gene_act=cell_gene.obs['gene_act_leiden_res0_3']
df_atac=cell_gene.obs['gene_act_leiden_res0_3']

df = pd.DataFrame({
    'gene_act': df_gene_act,
    'atac': df_atac
})
df['cell'] = df.index
gene_act_binary_matrix = pd.get_dummies(df.set_index('cell')['gene_act']).T
atac_binary_matrix = pd.get_dummies(df.set_index('cell')['atac']).T

def jaccard_similarity(x, y):
    intersection = np.logical_and(x, y).sum()
    union = np.logical_or(x, y).sum()
    return intersection / union if union != 0 else 0

num_cluster_gene_act = gene_act_binary_matrix.shape[0]
num_cluster_atac = atac_binary_matrix.shape[0]

jaccard_matrix = np.zeros((num_cluster_gene_act, num_cluster_atac))

for i in range(num_cluster_gene_act):
    for j in range(num_cluster_atac):
        jaccard_matrix[i, j] = jaccard_similarity(gene_act_binary_matrix.iloc[i], atac_binary_matrix.iloc[j])   

plt.figure(figsize=(10, 8))
sns.heatmap(jaccard_matrix, annot=True, cmap='coolwarm', 
            xticklabels=atac_binary_matrix.index, yticklabels=gene_act_binary_matrix.index, center=0)
plt.title('Jaccard Similarity between gene activity Clusters')
plt.xlabel('Clusters from gene activity clusters')
plt.ylabel('Clusters from gene activity clusters') 
plt.savefig(output_filename)    


end_time_read = time.time()
print("Time to process: "+str(end_time_read - start_time)+" seconds")



