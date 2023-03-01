# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from SEACells.core import summarize_by_SEACell
import anndata
from scipy.sparse import csr_matrix

# Some plotting aesthetics - will have to change as not in notebook
%matplotlib inline
sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100
def log_transform(ad, ps=0.1):
    ad.X.data = np.log2(ad.X.data + ps) - np.log2(ps)

# setting paths
os.chdir("/Users/hamrude/dev/repos/atac_neural_plate_border/local_test_data")
plot_path = "./plots/"
if not os.path.exists(plot_path):
    os.mkdir(plot_path)
    
"""
Load in test data
"""

#Load Data
#We recommend the use of scanpy Anndata objects as the preferred mode of loading and filtering data.

#RNA
# !mkdir data/
# !wget https://dp-lab-data-public.s3.amazonaws.com/SEACells/cd34_multiome_rna_with_labels.h5ad -O data/cd34_multiome_rna_with_labels.h5ad
# !wget https://dp-lab-data-public.s3.amazonaws.com/SEACells/cd34_multiome_atac_with_labels.h5ad -O data/cd34_multiome_atac_with_labels.h5ad
# Load the data using scanpy

"""

RNA Data Import and Processing -> will need to do this in R first

"""

# This object contains raw counts and the results of SEACells
#rna_ad = sc.read(f'data/cd34_multiome_rna_with_labels.h5ad')
rna_ad = sc.read("./AnnData_obj_RNA/rds_files/AnnData_metacells_assigned.h5ad") # ss8 test local

# check rna object by running UMAP
rna_ad

plt.figure(figsize=(8,8))
sc.pl.scatter(rna_ad, basis='umap', color="seurat_clusters", frameon=False)
plt.savefig(os.path.join(plot_path, "RNA_UMAP_clusters.png"))
plt.close()

# summarise by seacells
rna_meta_ad = summarize_by_SEACell(rna_ad, SEACells_label='SEACell')

### Processing ####
rna_meta_ad.obs['scHelper_cell_type'] = rna_ad.obs.groupby('SEACell').apply(lambda x: pd.Series(x['scHelper_cell_type']).mode())

# Normalize and log transform the data
sc.pp.normalize_per_cell(rna_meta_ad, counts_per_cell_after=10000)
log_transform(rna_meta_ad)

# Select highly variable genes - fails
sc.pp.highly_variable_genes(rna_meta_ad)

# Perform PCA and UMAPs
sc.tl.pca(rna_meta_ad, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(rna_meta_ad, n_neighbors=5)
sc.tl.umap(rna_meta_ad)
%matplotlib inline

## UMAP of metacells
sc.pl.scatter(rna_meta_ad, basis='umap', color='seurat_clusters', frameon=False)
sc.pl.scatter(rna_meta_ad, basis='umap', color='scHelper_cell_type', frameon=False)

## regressing?

"""

ATAC Data Import and Gene Scores Processing

"""

# We use the gene scores derived from ArchR aggregated by SEACell.

# read in data
#atac_ad = sc.read(f'data/cd34_multiome_atac_with_labels.h5ad')
atac_ad = sc.read("./SEACells_calculated/AnnData_metacells_assigned.h5ad")

# Check atac data by running UMAP - test data is ss8
atac_ad

plt.figure(figsize=(8,8))
sc.pl.scatter(atac_ad, basis='umap', color="clusters", frameon=False)
plt.savefig(os.path.join(plot_path, "ATAC_UMAP_clusters.png"))
plt.close()

atac_ad.obsm['GeneScores']
atac_ad.uns['GeneScoresColums']

# Extract gene scores matrix
genescores = pd.DataFrame(atac_ad.obsm['GeneScores'])
# why would the genes in the atac and rna neccessarily match like this??
#genescores.columns = rna_meta_ad.var_names
genescores.columns = atac_ad.uns['GeneScoresColums']
genescores.index = atac_ad.obs_names

# Aggregate gene scores by metacells
genescores = genescores.join(atac_ad.obs['SEACell']).groupby('SEACell').sum()

# Create Anndata
atac_meta_ad = anndata.AnnData(genescores)
atac_meta_ad.obs = atac_ad.obs.loc[atac_meta_ad.obs_names]
atac_meta_ad.X = csr_matrix(atac_meta_ad.X)

# Normalize the ATAC counts analogous to RNA counts
sc.pp.normalize_per_cell(atac_meta_ad, counts_per_cell_after=10000)
log_transform(atac_meta_ad)

#############   REACHED HERE
# Compute PCA using the highly variable genes from RNA 
hvg = atac_meta_ad.var_names.intersection(rna_meta_ad.var_names[rna_meta_ad.var['highly_variable']])
atac_hvg = pd.Series(False, index=atac_meta_ad.var_names)
atac_hvg[hvg] = True
atac_meta_ad.var['highly_variable'] = atac_hvg
sc.tl.pca(atac_meta_ad, n_comps=50)

# ATAC gene scores umap
sc.pp.neighbors(atac_meta_ad, use_rep='X_pca', n_neighbors=5)
sc.tl.umap(atac_meta_ad)
sc.pl.scatter(atac_meta_ad, basis='umap', color='celltype', frameon=False)

"""

Integrating ATAC and RNA

"""

# Annotate the Anndata objects based on their dataset
atac_meta_ad.obs['Dataset'] = 'ATAC'
rna_meta_ad.obs['Dataset'] = 'RNA'

# Create a combined Anndata object only including hvg
hvg = atac_meta_ad.var_names[atac_meta_ad.var['highly_variable']].intersection(rna_meta_ad.var_names)
comb_ad = rna_meta_ad[:, hvg].concatenate(atac_meta_ad[:, hvg])
/usr/local/anaconda3/envs/seacells/lib/python3.8/site-packages/anndata/_core/anndata.py:1785: FutureWarning: X.dtype being converted to np.float32 from float64. In the next version of anndata (0.9) conversion will not be automatic. Pass dtype explicitly to avoid this warning. Pass `AnnData(X, dtype=X.dtype, ...)` to get the future behavour.
  [AnnData(sparse.csr_matrix(a.shape), obs=a.obs) for a in all_adatas],
# PCA 
from sklearn.decomposition import TruncatedSVD
svd = TruncatedSVD(n_components=30)
comb_ad.obsm['X_pca'] = svd.fit_transform(comb_ad.X)
Integration using linear transport
Inputs are

xs, Matrix of "source" RNA expression levels from scRNA-seq
xt, Matrix of "target" gene activities calcualated from ATAC-seq peak accessibility
The number of metacells can be variable between the two modalities but the number of columns should remain the same. We will use the jointly computed PCs as input for integration

Optional paramters

ws: Vector of weights representing the size of each cluster/metacell from RNA
wt: Same as above but for ATAC
rho: Float in range [0,1] representing whether the final transformation should be closer to RNA distribution (0) or ATAC distribution (1). Default value is 1.
reg: Small float to make sure covariance matrices are invertible
xs = comb_ad[comb_ad.obs['Dataset'] == 'RNA'].obsm['X_pca']
xt = comb_ad[comb_ad.obs['Dataset'] == 'ATAC'].obsm['X_pca']
from SEACells.domainadapt import LinearOT

model = LinearOT(rho=0)

xs_transformed, xt_transformed = model.fit_transform(xs, xt)
comb_ad.obsm['X_pca_transformed'] = pd.DataFrame(xs_transformed).append(pd.DataFrame(xt_transformed)).values
/var/folders/x3/p87l8d1n5qvcxb8hjggq4cg80000gq/T/ipykernel_84957/2157489645.py:1: FutureWarning: The frame.append method is deprecated and will be removed from pandas in a future version. Use pandas.concat instead.
  comb_ad.obsm['X_pca_transformed'] = pd.DataFrame(xs_transformed).append(pd.DataFrame(xt_transformed)).values
Integration result
sc.pp.neighbors(comb_ad, use_rep='X_pca_transformed', n_neighbors=5)
sc.tl.umap(comb_ad)
sc.pl.scatter(comb_ad, basis='umap', color=['celltype', 'Dataset'], frameon=False)

 
Mapping ATAC SEACells to RNA SEACellss
# Mapping ATAC SEACells to RNA SEACellss
rna_PCs = pd.DataFrame(xs_transformed)
atac_PCs = pd.DataFrame(xt_transformed)

rna_PCs.index = rna_meta_ad.obs_names
atac_PCs.index = atac_meta_ad.obs_names
def compute_distances(df1, df2):
    result = pd.DataFrame(np.zeros((df1.shape[0],df2.shape[0]))).set_index(df1.index)
    result.columns = df2.index

    for mc1, row1 in df1.iterrows():
        for mc2, row2 in df2.iterrows():
            result.loc[mc1, mc2] = np.sqrt(((row1-row2)**2).sum())

    return result

def get_match_ranks(t,k=5):
    """
    Compute mutually ranked metacell matches.
    """
    arr = t.values.T
    ind = np.argsort(arr, axis=0)[:k].T

    df = pd.DataFrame(t.columns.values[ind])
    df.index = t.index
    
    return df

def MNN(rna_PCs, atac_PCs, k=3):
    """
    Construct a mapping between RNA and ATAC SEACells, if an RNA SEACell and 
    ATAC SEACell are mutually in each others' k-nearest neighbors.
    
    Nearest neighbors are computed based on distance in transformed PC space,
    following domain adaption by optimal transport.
    """
    
    # Compute distances between PCs in each data set
    dists = compute_distances(rna_PCs, atac_PCs)
    
    # Find the top k nearest neighbors in each row of the distance matrix
    r1 = get_match_ranks(dists, k)
    r2 = get_match_ranks(dists.T, k)
    # Add edges between samples if they appear in each others top-ranks 

    edges = []
    for index, row in r1.iterrows():
        for nbr in row.values:
            if index in r2.loc[nbr].values:
                edges.append((index, nbr))

    edges = pd.DataFrame(edges)
    edges.columns = ['RNA', 'ATAC']
    return edges
mapping = MNN(rna_PCs, atac_PCs)
mapping.head()
RNA	ATAC
0	cd34_multiome_rep2#CATTCCTCACCCACAG-1	cd34_multiome_rep1#GGCGGTAAGTAGCGGG-1
1	cd34_multiome_rep2#CATTCCTCACCCACAG-1	cd34_multiome_rep1#TACAGCTAGACCATAC-1
2	cd34_multiome_rep2#CAATCGCCACGTAAGG-1	cd34_multiome_rep2#TAACCTAAGGATTGAG-1
3	cd34_multiome_rep1#CCAAGTTAGGACCTCA-1	cd34_multiome_rep1#ATGGCTGTCAATGTCA-1
4	cd34_multiome_rep1#CCAAGTTAGGACCTCA-1	cd34_multiome_rep1#ATTGCACAGTCAATTG-1
 