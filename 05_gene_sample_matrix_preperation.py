import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


adata=sc.read('../result/1_preprocess/adata/discovery_celltypist_annotated.h5ad')
adata.obs = adata.obs.rename(columns={'predicted_labels':'celltype'})
raw=sc.read('../result/1_preprocess/adata/inner_QC.h5ad')
raw.obs=adata.obs

del(adata)

#Discovery cohort
adata=raw
for ct in adata.obs['celltype'].unique():
    temp = adata[adata.obs['celltype'] == ct].copy()
    temp.obs['sample'] = temp.obs['sample'].astype('category')
    res = pd.DataFrame(columns=temp.var_names, index=temp.obs['sample'].cat.categories)
    for clust in temp.obs['sample'].cat.categories: 
        res.loc[clust] = temp[temp.obs['sample'].isin([clust]),:].X.mean(0)
    res.to_csv('../result/1_preprocess/pseudobulk/discovery_cohort/{}.csv'.format(ct))


#Validation cohort
adata  = sc.read('../result/1_preprocess/adata/validation/validation_celltypist_annotated.h5ad')
adata.obs = adata.obs.rename(columns={'predicted_labels':'celltype'})

for ct in adata.obs['celltype'].unique():
    temp = adata[adata.obs['celltype'] == ct].copy()
    temp.obs['sample'] = temp.obs['sample'].astype('category')
    res = pd.DataFrame(columns=temp.var_names, index=temp.obs['sample'].cat.categories)
    for clust in temp.obs['sample'].cat.categories: 
        res.loc[clust] = temp[temp.obs['sample'].isin([clust]),:].X.mean(0)
    res.to_csv('../result/1_preprocess/pseudobulk/validation_cohort/{}.csv'.format(ct))