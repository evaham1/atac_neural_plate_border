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

# setting paths
os.chdir("/Users/hamrude/dev/repos/atac_neural_plate_border/local_test_data")
plot_path = "./seacell_integration/plots/"
if not os.path.exists(plot_path):
    os.mkdir(plot_path)
    
"""
Load in data
"""

rna_ad = sc.read("./test_inputs/test_input_AnnData_seacells_processed_for_integration/AnnData_ss8_RNA.h5ad") # ss8 test local

# check rna object by running UMAP
rna_ad

#plt.figure(figsize=(8,8))
sc.pl.scatter(rna_ad, basis='umap', color="scHelper_cell_type", frameon=False)

rna_ad.obs['cell_colours']

sc.pl.umap(adata, color='cell_type', color_map=cell_type_colors)

sc.pl.scatter(rna_ad, basis='umap', color="scHelper_cell_type", frameon=False)
#plt.savefig(os.path.join(plot_path, "RNA_UMAP_clusters.png"))
#plt.close()

atac_ad = sc.read("./test_inputs/test_input_AnnData_seacells_processed_for_integration/AnnData_ss8_ATAC.h5ad") # ss8 test local

# check rna object by running UMAP
rna_ad

#plt.figure(figsize=(8,8))
sc.pl.scatter(atac_ad, basis='umap', color="scHelper_cell_type", frameon=False)
#plt.savefig(os.path.join(plot_path, "ATAC_UMAP_clusters.png"))
#plt.close()


"""

Preparing combined AnnData object

"""

# Annotate the Anndata objects based on their dataset
atac_ad.obs['Dataset'] = 'ATAC'
rna_ad.obs['Dataset'] = 'RNA'

# Identify highly variable genes in RNA object
sc.pp.highly_variable_genes(rna_ad)

# Identify list of hvg which are also present in atac object
hvg = rna_ad.var_names[rna_ad.var['highly_variable']].intersection(atac_ad.var_names)
hvg

# Create a combined Anndata object only including hvg
comb_ad = rna_ad[:, hvg].concatenate(atac_ad[:, hvg])

comb_ad

# PCA 
from sklearn.decomposition import TruncatedSVD
svd = TruncatedSVD(n_components=30)
comb_ad.obsm['X_pca'] = svd.fit_transform(comb_ad.X)


"""

Running integration

"""

xs = comb_ad[comb_ad.obs['Dataset'] == 'RNA'].obsm['X_pca']
xt = comb_ad[comb_ad.obs['Dataset'] == 'ATAC'].obsm['X_pca']


from SEACells.domainadapt import LinearOT

model = LinearOT(rho=0)

xs_transformed, xt_transformed = model.fit_transform(xs, xt)

comb_ad.obsm['X_pca_transformed'] = pd.DataFrame(xs_transformed).append(pd.DataFrame(xt_transformed)).values

"""

View integration results

"""

sc.pp.neighbors(comb_ad, use_rep='X_pca_transformed', n_neighbors=5)

sc.tl.umap(comb_ad)

sc.pl.scatter(comb_ad, basis='umap', color=['scHelper_cell_type', 'Dataset'], frameon=False)






