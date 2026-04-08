import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import doubletdetection
import os
import scipy.sparse as sp
from scipy.sparse import csr_matrix


"""
Preprocess scRNA seq raw data for 28 different data source. The preprocessing steps include:
1. Renaming files to match the format required by `sc.read_10x_mtx()` and loading the data.
2. Assigning a unique identifier to each sample and adding age and sex information.
3. Performing doublet detection and saving the results.
4. Merging all samples and saving the combined dataset as `raw.h5ad`. The dataset includes metadata such as age, sex, pub, sample, doublet, and doublet_score, while the matrix content remains unchanged.

The following example demonstrates the process for sample P23.
"""

# Initialize the doublet detection classifier
clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)

# Define a function to rename raw files: `rename_raw()`
def rename_raw(directory):
    """
    Rename files in the directory to match the standard format required by `sc.read_10x_mtx()`.
    The function ensures the presence of `matrix.mtx`, `barcodes.tsv`, and `features.tsv`.

    Parameters:
        directory (str): Path to the directory containing the raw files.
    """
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        
        if 'matrix.mtx' in filename and not filename == 'matrix.mtx':
            os.rename(filepath, os.path.join(directory, 'matrix.mtx.gz'))
            print(f"Renamed {filename} to matrix.mtx")
        
        elif 'barcodes.tsv' in filename and not filename == 'barcodes.tsv':
            os.rename(filepath, os.path.join(directory, 'barcodes.tsv.gz'))
            print(f"Renamed {filename} to barcodes.tsv")
        
        elif 'features.tsv' in filename and not filename == 'features.tsv':
            os.rename(filepath, os.path.join(directory, 'features.tsv.gz'))
            print(f"Renamed {filename} to features.tsv")
            
        elif 'genes.tsv' in filename and not filename == 'genes.tsv':
            os.rename(filepath, os.path.join(directory, 'genes.tsv.gz'))
            print(f"Renamed {filename} to features.tsv")

# Age and sex information (retrieved from the publication)
age_dict = {'HC1': 90, 'HC2': 68, 'HC3': 38, 'HC4': 84, 'HC5': 70}
sex_dict = {'HC1': 'Male', 'HC2': 'Female', 'HC3': 'Male', 'HC4': 'Female', 'HC5': 'Female'}

print("Starting the loop for processing samples...")

# Process each sample
for pub in ['P23']:
    for sample in ['HC1', 'HC2', 'HC3', 'HC4', 'HC5']:
        print(f"Processing sample: {sample}   ", end="")
        
        # Rename files to match the required format
        rename_raw(f"../data/adata/rawdata/P23/{sample}")
        
        # Load the dataset
        adata = sc.read_10x_mtx(f"../data/adata/raw_data/P23/{sample}/")
        
        # Add metadata to the dataset
        adata.obs['pub'] = pub  # Add publication identifier (used as batch information)
        adata.obs['sample'] = f"{pub}_{sample}"  # Add sample identifier
        adata.obs_names = [f"{pub}_{sample}_{name}" for name in adata.obs_names]  # Rename cells to avoid duplication
        adata.obs['age'] = age_dict[sample]  # Add age information
        adata.obs['sex'] = sex_dict[sample]  # Add sex information
        
        # Perform doublet detection
        doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
        doublet_score = clf.doublet_score()
        adata.obs["doublet"] = doublets  # Save doublet detection results
        adata.obs["doublet_score"] = doublet_score
        
        # Merge datasets
        if sample == 'HC1':
            adata_merge = adata
        else:
            adata_merge = ad.concat([adata_merge, adata], join='outer')  # Use outer join to include all genes

# Save the merged dataset
adata_merge.write('../data/adata/raw_data/P23_raw.h5ad')


