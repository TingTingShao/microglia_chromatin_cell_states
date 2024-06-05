#!/usr/bin/env python
import pandas as pd

bed_file_path = '/lustre1/project/stg_00079/students/tingting/data/GWAS/020results_expansion/Bellenguez_etal_Stage1_result.expanded.bed'
df=pd.read_csv(bed_file_path, sep="\t", header=None)
df.columns=['Chromosome', 'start_pos', 'snp_pos', 'rsid', 'r_square', 
            'Chromosome_AD', 'start_pos_AD', 'snp_pos_AD', 'rsid_AD', 'p_value']

df.to_parquet("/lustre1/project/stg_00079/students/tingting/data/GWAS/020results_expansion/Bellenguez_etal_Stage1_result.expanded.parquet", index=False)