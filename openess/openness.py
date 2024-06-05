#!/usr/bin/env python
import anndata
# import snapatac2 as snap
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix, issparse

adat=anndata.read_h5ad('/data/leuven/351/vsc35107/lustre1_stt/data/sun/snap2_allfragments/08/microglia_1.h5ad')

var_names=adat.var_names.tolist()
X_dense = adat.X

split_data = []
for item in var_names:
    chr_part, range_part = item.split(':')
    start, end = range_part.split('-')
    split_data.append([chr_part, int(start), int(end)])

def calculate_column_averages(matrix):
    if issparse(matrix):
        return matrix.mean(axis=0).A1  # .A1 to convert to 1D array
    else:
        return np.mean(matrix, axis=0).flatten()

column_averages_dense = calculate_column_averages(X_dense)

df = pd.DataFrame(split_data, columns=['chr', 'start', 'end'])
df['openness'] = column_averages_dense

df.to_parquet('/data/leuven/351/vsc35107/lustre1_stt/data/sun/snap2_allfragments/08/openness.parquet')



