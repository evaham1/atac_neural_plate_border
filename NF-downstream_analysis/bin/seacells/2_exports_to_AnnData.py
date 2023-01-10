#!/usr/bin/env python

# # load libraries
import os
import pandas as pd
import numpy as np

import scanpy as sc
import pyranges as pr
import warnings
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import SEACells

print("Libraries loaded")

# Choose which GUI matplotlib uses to print out plots - for notebook only
#matplotlib.use('TkAgg') # prints them in a separate window
#%matplotlib notebook

# set data path and check contents
data_dir = os.path.expanduser('./ArchR_exported_data/')
arr = os.listdir(data_dir)
print(arr)

# Counts data - sparse COO matrix
from scipy.io import mmread
counts = mmread(data_dir + 'peak_counts/counts.mtx')
print(counts.todense()[:10])

# Cell information
cells = pd.read_csv(data_dir + 'peak_counts/cells.csv', index_col=0).iloc[:, 0]
print(cells.head())

# Peaks information
peaks = pd.read_csv(data_dir + 'peak_counts/peaks.csv', index_col=0)
peaks.index = peaks['seqnames'] + ':' + peaks['start'].astype(str) + '-' + peaks['end'].astype(str)
print(peaks.head())

# Make AnnData object
ad = sc.AnnData(counts.T)
ad.obs_names = cells
ad.var_names = peaks.index
for col in peaks.columns:
    ad.var[col] = peaks[col]
ad.X = ad.X.tocsr()
print(ad)

# Read in reduced dimensions and add to AnnData
ad.obsm['X_svd'] = pd.read_csv(data_dir + 'svd.csv', index_col=0).loc[ad.obs_names, : ].values
print(ad.obsm)

# Read in cell metadata
cell_meta = pd.read_csv(data_dir + 'cell_metadata.csv', index_col=0).loc[ad.obs_names, : ]
for col in cell_meta.columns:
    ad.obs[col] = cell_meta[col].values
print(ad.obs)

# Gene scores
gene_scores = pd.read_csv(data_dir + 'gene_scores.csv', index_col=0).T

ad.obsm['GeneScores'] = gene_scores.loc[ad.obs_names, :].values
ad.uns['GeneScoresColums'] = gene_scores.columns.values

# Leiden and UMAP
warnings.filterwarnings('ignore')
sc.pp.neighbors(ad, use_rep='X_svd')
sc.tl.umap(ad)
sc.tl.leiden(ad)
warnings.filterwarnings('default')

# Plots to check - need to get these to print??
sc.pl.scatter(ad, basis='umap', color= "#fe57a1")

sc.pl.scatter(ad, basis='umap', color='leiden')

# Save AnnData object
ad.write(data_dir + '../ss8.h5ad')