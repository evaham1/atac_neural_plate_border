#!/usr/bin/env python

# load libraries
import os
import argparse
import sys
import pandas as pd
import numpy as np

import scanpy as sc
import pyranges as pr
import warnings
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import SEACells
from scipy.io import mmread

def parse_args(args=None):
    Description = " "
    Epilog = " "

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help="Input path.", metavar='')
    parser.add_argument('-c', '--ncores', type=str, help="Number of CPUs.", metavar='')
    return parser.parse_args(args)

def read_input(input):
    counts = mmread(input + '/peak_counts/counts.mtx')
    print("Count data loaded!")
    cells = pd.read_csv(input + '/peak_counts/cells.csv', index_col=0).iloc[:, 0]
    print("Cell data loaded!")
    peaks = pd.read_csv(input + '/peak_counts/peaks.csv', index_col=0)
    print("Peak data loaded!")
    return(counts, cells, peaks)

def make_Anndata(counts, cells, peaks):
    ad = sc.AnnData(counts.T)
    ad.obs_names = cells
    ad.var_names = peaks.index
    for col in peaks.columns:
        ad.var[col] = peaks[col]
    ad.X = ad.X.tocsr()
    return(ad)

def cluster_and_UMAP(ad):
    warnings.filterwarnings('ignore')
    sc.pp.neighbors(ad, use_rep='X_svd')
    print("Neighbours found!")
    sc.tl.umap(ad)
    print("UMAP calculated!")
    sc.tl.leiden(ad)
    print("Clustered!")
    warnings.filterwarnings('default')
    return(ad)

def main(args=None):

    # read in command line arguments
    args = parse_args(args)

    # set output paths
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)
    if not os.path.exists(rds_path):
        os.mkdir(rds_path)

    # Load data
    counts, cells, peaks = read_input(args.input)

    # add peak indices
    peaks.index = peaks['seqnames'] + ':' + peaks['start'].astype(str) + '-' + peaks['end'].astype(str)
    #print(peaks.head())

    # Make AnnData object
    ad = make_Anndata(counts, cells, peaks)
    print("AnnData object created!")
    #print(ad)

    # Read in reduced dimensions and add to AnnData
    ad.obsm['X_svd'] = pd.read_csv(args.input + '/svd.csv', index_col=0).loc[ad.obs_names, : ].values
    print("Reduced dims added!")
    #print(ad.obsm)

    # Read in cell metadata
    cell_meta = pd.read_csv(args.input + '/cell_metadata.csv', index_col=0).loc[ad.obs_names, : ]
    for col in cell_meta.columns:
        ad.obs[col] = cell_meta[col].values
    print("Cell metadata added!")
    #print(ad.obs)

    # Read in gene scores
    gene_scores = pd.read_csv(args.input + '/gene_scores.csv', index_col=0).T
    ad.obsm['GeneScores'] = gene_scores.loc[ad.obs_names, :].values
    ad.uns['GeneScoresColums'] = gene_scores.columns.values
    print("Gene scores added!")

    # Custering and UMAP calculate
    ad = cluster_and_UMAP(ad)

    # Print UMAP plot to check data looks like before
    fig = sc.pl.scatter(ad, basis='umap', color= 'leiden')
    plt.savefig(os.path.join(plot_path, "UMAP.png"))

    # Save AnnData object
    ad.write(os.path.join(rds_path, 'AnnData.h5ad'))

    # read in anndata for quicker debugging
    #ad = sc.read('AnnData.h5ad')
    #print(ad)

if __name__ == '__main__':
    sys.exit(main())