import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import doubletdetection
import os
import scanpy.external as sce

# Integrate multiple files
P_num = {'P0', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15',
         'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'P24', 'P25', 'P26', 'P27', 'P28'}

for i in P_num:
    temp = sc.read_h5ad('../data/adata/raw_data/{}_raw.h5ad'.format(i))
    
    temp.obs['age'] = temp.obs['age'].cat.codes.astype(float)
    adata_inner.obs['age'] = adata_inner.obs['age'].astype('category')
    adata_inner.obs['age'] = adata_inner.obs['age'].cat.codes.astype(float)
    print(adata_inner.obs['age'].dtype)
    
    # Concatenate datasets, keeping only shared genes across all 28 samples
    adata_inner = ad.concat([adata_inner, temp], join='inner')
    print(i)

adata_inner.write('../result/1_preprocess/raw_data/inner_raw.h5ad')

# Load the integrated dataset into memory
adata = sc.read('../result/1_preprocess/raw_data/inner_raw.h5ad', backed='r').to_memory()
adata.X = adata.X.astype(np.float32)

# Quality control
adata = adata[adata.obs['doublet'] == 0]
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Annotate mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 25, :]

adata.write('../result/1_preprocess/preprocessed/inner_QC.h5ad')

# Normalization
sc.pp.normalize_total(adata, target_sum=1e6)
adata.write('../result/1_preprocess/preprocessed/inner_QC_norm.h5ad')

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="pub")
adata = adata[:, adata.var.highly_variable]

# Scale the data
sc.pp.scale(adata, max_value=10)

# Dimensionality reduction
sc.tl.pca(adata, svd_solver='arpack')

# Batch correction using Harmony
sce.pp.harmony_integrate(adata, key=['pub'])  # Harmony is used to remove batch effects
print('Harmony integration is complete!')

# Clustering 
n_n = 35
n_pc = 27

sc.pp.neighbors(adata, n_neighbors=n_n, n_pcs=n_pc, use_rep='X_pca_harmony')
sc.tl.umap(adata)

adata.write('../result/1_preprocess/preprocessed/harmony_inner.h5ad')