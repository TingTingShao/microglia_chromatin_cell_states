#!/usr/bin/env python

"""
visualise the cellular proportions before and after batch correction
03.23 sample correction (harmony vs. MNC), cellular proportion in each cluster (resolution 0.3) vs. sample, region, subject, ad2

run example: cell_proportion.py -a microglia_1.h5ads -r 0.5 -o 00.png 01.png 02.png 03.png
"""

import snapatac2 as snap
import os, sys, argparse
import matplotlib.pyplot as plt
import numpy as np
import duckdb as db
import pandas as pd
import time
import seaborn as sns

start_time = time.time()

input_filename = None
output_filenames = None
resolution=0.3

parser=argparse.ArgumentParser(description="visualise the cellular proportion", \
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-a", '--input_file', nargs=1, default=None, help="input the anndataset")
parser.add_argument("-r", '--resolution', nargs=1, default=0.3, type=float, help="input the anndataset")
parser.add_argument("-o", '--output_files', nargs='+', default=None, help="output file names")

args=parser.parse_args()
print(args)

input_filename = args.input_file[0]
# print(input_filename)
resolution=args.resolution[0]
output_filenames = args.output_files
# print(resolution)
print("Output files:", output_filenames)

dat=snap.read_dataset(input_filename)

db.register('datobs', dat.obs[:])

def generate_combined_heatmap(category, output_filename):
    column_prefix = "leiden_mnc"
    # Construct the query using the determined column prefix
    query = f"""
    SELECT
        {category},
        "{column_prefix}_{resolution}",
        COUNT(*) AS nuclei_count
    FROM
        datobs
    GROUP BY
        {category},
        "{column_prefix}_{resolution}"
    ORDER BY
        {category},
        "{column_prefix}_{resolution}";
    """
    df = db.query(query).df()
    # heatmap_data = df.pivot_table(index=f"{column_prefix}_{resolution}", columns=category, values='nuclei_count', fill_value=0)
    # plt.figure(figsize=(20, 10))
    # sns.heatmap(heatmap_data, cmap="YlGnBu", linewidths=.5)
    # plt.title(f'{category} Count Distribution in Each Cluster for cluster')
    # plt.xlabel(f'{category.capitalize()}')
    # plt.ylabel('Cluster')
    # plt.savefig(output_filename)
    # plt.close()

    plt.figure(figsize=(12, 8))
    sns.barplot(data=df, x=f"{column_prefix}_{resolution}", y='nuclei_count', hue=category, estimator=sum, ci=None)
    plt.title(f'{category} Distribution of cellular counts in each cluster')
    plt.xlabel('Cluster')
    plt.ylabel('Nuclei Count')
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.savefig(output_filename)
    plt.close()


categories = ['sample', 'region', 'subject', 'ad']
# pre_correction_filenames = [f"pre_{filename}" for filename in output_filenames]
# post_correction_filenames = [f"post_{filename}" for filename in output_filenames]

# Loop through categories and generate heatmaps for both pre- and post-correction
if len(output_filenames) != len(categories):
    print("Error: The number of output filenames does not match the number of categories.")
else:
    for category, post_filename in zip(categories, output_filenames):
        # generate_combined_heatmap(category, pre_filename, correction_state='pre')  # For pre-correction
        generate_combined_heatmap(category, post_filename)  # For post-correction


dat.close()

end_time_read = time.time()
print("Time to process: "+str(end_time_read - start_time)+" seconds")
