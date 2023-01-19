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

# Some plotting aesthetics
sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100

def parse_args(args=None):
    Description = " "
    Epilog = " "

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help="Input path.", metavar='')
    parser.add_argument('-c', '--ncores', type=str, help="Number of CPUs.", metavar='')
    return parser.parse_args(args)

def main(args=None):

    # read in command line arguments
    args = parse_args(args)
    print(args.input)

    # set output paths
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)
    if not os.path.exists(rds_path):
        os.mkdir(rds_path)

    # Load data
    ad = sc.read(args.input + '/AnnData.h5ad')
    print(ad)

    # Plot cell-types for reference
    plt.figure(figsize=(8,8))
    sc.pl.scatter(ad, basis='umap', color='stage_scHelper_cell_type_old', frameon=False)
    plt.savefig(os.path.join(plot_path, "UMAP_stage_scHelper_cell_type_old.png"))
    plt.close()

    plt.figure(figsize=(8,8))
    sc.pl.scatter(ad, basis='umap', color='clusters', frameon=False)
    plt.savefig(os.path.join(plot_path, "UMAP_clusters.png"))
    plt.close()

    #####   Run SEACells
    plot_path = "./plots/computing_SEACells/"
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)

    # Set parameters 
    #n_SEACells = 10 # for s8 data that has 7,409 cells should use 90, use less for now to speed up
    n_SEACells = 100 # for TL data that has 80,232 cells should be 1000 to get one metacell per 80 cells
    build_kernel_on = 'X_svd' 
    n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells

    # Build SEACells model
    print("initializing SEACElls model...")
    model = SEACells.core.SEACells(ad, 
                  build_kernel_on=build_kernel_on, 
                  n_SEACells=n_SEACells, 
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5)
    model.construct_kernel_matrix()
    M = model.kernel_matrix
    print("SEACElls model made!")

    # Plot correlation heatmap
    plt.figure(figsize=(8,8))
    sns.clustermap(M.toarray()[:500,:500])
    plt.savefig(os.path.join(plot_path, "Correlation_heatmap.png"))
    plt.close()

    # Initialize archetypes
    print("Initializing archetypes...")
    model.initialize_archetypes()
    print("Archetypes initialized!")

    # Plot the initilization to ensure they are spread across phenotypic space
    SEACells.plot.plot_initialization(ad, model, save_as = os.path.join(plot_path, "Initialisations_UMAP.png"))

    # Iterate with different initialisations
    print("Fitting SEACElls model...")
    model.fit(min_iter=10, max_iter=100) # for TL data
    #model.fit(min_iter=1, max_iter=20) # to speed up for testing
    print(f'Ran for {len(model.RSS_iters)} iterations')
    for _ in range(5):
        model.step()
    print(f'Ran for {len(model.RSS_iters)} iterations')
    print("SEACells model fitted!")

    # Check for convergence 
    model.plot_convergence(save_as = os.path.join(plot_path, "Convergence_plot.png"))

    ## Visualise cell assignments
    plot_path = "./plots/soft_SEACell_assignments/"
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)

    # Visualise hard assignments
    ad.obs[['SEACell']].head()
    model.get_hard_assignments().head()

    # Visualise soft assignments
    plt.figure(figsize=(6,4))
    sns.distplot((model.A_.T > 0.1).sum(axis=1), kde=False)
    plt.title(f'Non-trivial (> 0.1) assignments per cell')
    plt.xlabel('# Non-trivial SEACell Assignments')
    plt.ylabel('# Cells')
    plt.savefig(os.path.join(plot_path, "Non-trivial_SEACell_Assignments.png"), bbox_inches = 'tight')
    plt.close()

    plt.figure(figsize=(10,10)) # creates new figure
    b = np.partition(model.A_.T, -5)    
    sns.heatmap(np.sort(b[:,-5:])[:, ::-1], cmap='viridis', vmin=0)
    plt.title('Strength of top 5 strongest assignments')
    plt.xlabel('$n^{th}$ strongest assignment')
    plt.savefig(os.path.join(plot_path, "Strength_of_top_5_strongest_assignments.png"), dpi=199)
    plt.close()

    labels,weights = model.get_soft_assignments()

    #####   Summarise counts by metacells
    # summarising needs to be run on a layer so copy the matrix to a layer called raw
    ad.layers["raw"] = ad.X
    ad.layers.keys()

    # summarise counts by hard labels -> SEACell_ad
    print("Summarising counts for hard labels...")
    SEACell_ad = SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell', summarize_layer='raw')
    SEACell_ad.obs.head()

    # summarise counts by soft labels -> SEACell_soft_ad
    # print("Summarising counts for soft labels...")
    # SEACell_soft_ad = SEACells.core.summarize_by_soft_SEACell(ad, model.A_, celltype_label='stage_scHelper_cell_type_old',summarize_layer='raw', minimum_weight=0.05)
    # SEACell_soft_ad.obs.head()

    ### Do we need to normalise??
    # sc.pp.normalize_total
    # sc.pp.log1p

    #####    Visualising metacells
    plot_path = "./plots/visualise_SEACells/"
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)

    # where metacells are on UMAP
    SEACells.plot.plot_2D(ad, key='X_umap', colour_metacells=False, save_as = os.path.join(plot_path, "Metacells_on_UMAP.png"))
    SEACells.plot.plot_2D(ad, key='X_umap', colour_metacells=True, save_as = os.path.join(plot_path, "Metacells_on_UMAP_coloured.png"))

    # sizes of metacells
    SEACells.plot.plot_SEACell_sizes(ad, bins=5, save_as = os.path.join(plot_path, "SEACell_sizes.png"))

    #####    Checking quality of metacells
    plot_path = "./plots/QC_SEACells/"
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)

    # purity of metacells - clusters on full data
    SEACell_purity = SEACells.evaluate.compute_celltype_purity(ad, 'clusters')
    SEACell_purity.head()
    plt.figure(figsize=(5,5))
    sns.boxplot(data=SEACell_purity, y='clusters_purity')
    plt.title('ArchR clusters Purity')
    sns.despine()
    plt.savefig(os.path.join(plot_path, "Metacell_ArchR_clusters_purity.png"), bbox_inches="tight")
    plt.close()

    # purity of metacells - scHelper cell type from individual stages
    SEACell_purity = SEACells.evaluate.compute_celltype_purity(ad, 'stage_scHelper_cell_type_old')
    SEACell_purity.head()
    plt.figure(figsize=(5,5))
    sns.boxplot(data=SEACell_purity, y='stage_scHelper_cell_type_old_purity')
    plt.title('scHelper cell type Purity')
    sns.despine()
    plt.savefig(os.path.join(plot_path, "Metacell_stage_scHelper_cell_type_purity.png"), bbox_inches="tight")
    plt.close()

    # purity of metacells - clusters on individual stages
    SEACell_purity = SEACells.evaluate.compute_celltype_purity(ad, 'stage_clusters')
    SEACell_purity.head()
    plt.figure(figsize=(5,5))
    sns.boxplot(data=SEACell_purity, y='stage_clusters_purity')
    plt.title('ArchR clusters Purity')
    sns.despine()
    plt.savefig(os.path.join(plot_path, "Metacell_ArchR_stage_clusters_purity.png"), bbox_inches="tight")
    plt.close()
    
    # compactness
    compactness = SEACells.evaluate.compactness(ad, 'X_svd')
    plt.figure(figsize=(5,5))
    sns.boxplot(data=compactness, y='compactness')
    plt.title('Compactness')
    sns.despine()
    plt.savefig(os.path.join(plot_path, "Compactness.png"), bbox_inches="tight")
    plt.close()

    # separation
    separation = SEACells.evaluate.separation(ad, 'X_svd',nth_nbr=1)
    plt.figure(figsize=(5,5))
    sns.boxplot(data=separation, y='separation')
    plt.title('Separation')
    sns.despine()
    plt.savefig(os.path.join(plot_path, "Separation.png"), bbox_inches="tight")
    plt.close()

    ### Saving data

    # Save AnnData object with metacell assignments
    ad.write(os.path.join(rds_path, 'AnnData_metacells_assigned.h5ad'))

    # Save summarised AnnData object 
    SEACell_ad.write(os.path.join(rds_path, 'AnnData_summarised_by_metacells.h5ad'))


if __name__ == '__main__':
    sys.exit(main())