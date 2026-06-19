import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import scrublet as scr
import argparse
import os
import scvi
from scvi.model import SCVI
import torch

def detect_and_remove_doublets(adata, expected_doublet_rate=0.075):
    """
    Detect doublets using Scrublet and annotate the AnnData object.

    Note: this function only computes and stores doublet scores and
    predictions in `adata.obs`. Actual removal/filtering is handled in
    `quality_control` so the user can visualize results before filtering.
    """
    # Ensure dense matrix for Scrublet; avoid forcing dense if already small
    matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)

    scrub = scr.Scrublet(matrix, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
    )

    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets

    return adata

def visualize_doublet_results(adata):
    """
    Visualize doublet detection results.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    plt.figure(figsize=(12, 4))

    plt.subplot(1, 2, 1)
    sns.histplot(adata.obs['doublet_score'], bins=50, kde=True)
    plt.axvline(x=0.25, color='r', linestyle='--', label='Threshold')
    plt.title('Distribution of Doublet Scores')
    plt.xlabel('Doublet Score')
    plt.ylabel('Cell Count')
    plt.legend()

    plt.subplot(1, 2, 2)
    sns.scatterplot(
        x=adata.obs['n_genes_by_counts'],
        y=adata.obs['doublet_score'],
        hue=adata.obs.get('predicted_doublet', None),
        alpha=0.5
    )
    plt.title('Doublet Scores vs. Gene Counts')
    plt.xlabel('Gene Counts per Cell')
    plt.ylabel('Doublet Score')

    plt.tight_layout()
    plt.savefig('doublet_analysis.png', dpi=300)
    plt.close()

def quality_control(adata, min_genes=400, max_genes=2500, max_mito_ratio=0.1):
    """
    Perform quality control on single-cell data, including doublet detection.
    Parameters:
        adata: AnnData object
        min_genes: Minimum gene count threshold per cell
        max_genes: Maximum gene count threshold per cell
        max_mito_ratio: Maximum mitochondrial gene expression ratio
    """
    print("Starting quality control...")
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-') 
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    
    adata = detect_and_remove_doublets(adata)
    visualize_doublet_results(adata)

    adata.obs['low_quality'] = (
        (adata.obs.n_genes_by_counts < min_genes)
        | (adata.obs.n_genes_by_counts > max_genes)
        | (adata.obs.pct_counts_mt > max_mito_ratio)
    )

    sc.pp.filter_genes(adata, min_cells=3)
    
    n_cells_before = adata.n_obs

    # Combine filters: remove predicted doublets and low-quality cells
    if 'predicted_doublet' in adata.obs.columns:
        combined_filter = (~adata.obs['predicted_doublet']) & (~adata.obs['low_quality'])
    else:
        combined_filter = (~adata.obs['low_quality'])

    adata = adata[combined_filter, :].copy()
    n_cells_after = adata.n_obs

    print(f"Cells before QC: {n_cells_before}, Cells after QC: {n_cells_after}")
    print(f"Filtered out {n_cells_before - n_cells_after} cells (doublets + low-quality).")

    return adata

def normalize_and_select_features(adata, n_top_genes=2000):
    """
    Normalize data and select highly variable genes.
    """
    print("Starting normalization and highly variable gene selection...")
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    sc.pp.highly_variable_genes(
        adata, 
        n_top_genes=n_top_genes,
        flavor='seurat_v3',
        batch_key='batch' if 'batch' in adata.obs.columns else None
    )
    
    print(f"Selected {n_top_genes} highly variable genes.")
    return adata

def scale_and_pca(adata):
    """
    Scale data and perform PCA.
    """
    print("Starting data scaling and PCA...")
    
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(adata, svd_solver='arpack')
    
    n_pcs = min(50, len(adata.obs))
    print(f"Using {n_pcs} principal components.")
    
    return adata, n_pcs

def batch_correction_with_scvi(adata, batch_key='batch', n_latent=30):
    """
    Perform batch correction using scVI.
    """
    print("Starting batch correction (scVI)...")
    
    # Ensure data is in log space
    if 'log1p' not in adata.uns:
        sc.pp.log1p(adata)

    # Prepare scVI-compatible AnnData and train the model
    adata_scvi = adata.copy()
    scvi.data.setup_anndata(adata_scvi, batch_key=batch_key)

    model = SCVI(adata_scvi, n_latent=n_latent)
    model.train(max_epochs=100, use_gpu=torch.cuda.is_available())

    # Retrieve batch-corrected latent representation and store in obsm
    latent = model.get_latent_representation()
    adata.obsm['X_scVI'] = latent

    print("Batch correction completed.")
    return adata

def cluster_and_visualize(adata, n_neighbors=15, resolution=1.0):
    """
    Perform cell clustering and UMAP visualization.
    """
    print("Starting cell clustering and UMAP visualization...")
    
    sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors=n_neighbors)
    sc.tl.louvain(adata, resolution=resolution)
    sc.tl.umap(adata, min_dist=0.3)

    n_clusters = adata.obs['louvain'].nunique()
    print(f"Identified {n_clusters} cell clusters.")

    return adata


def main():
    """
    Main function: Execute the full analysis workflow and save intermediate results at key steps.
    """
    INPUT_PATH = "../result/1_preprocess/raw_data/PJ_raw.h5ad"
    QC_PATH = "../result/1_preprocess/validation/PJ_QC.h5ad"
    NORMALIZED_PATH = "../result/1_preprocess/validation/PJ_QC_NORM.h5ad"
    BATCH_CORRECTED_PATH = "../result/1_preprocess/validation/PJ_QC_NORM_batch.h5ad"
    CLUSTERED_PATH = "../result/1_preprocess/validation/PJ_QC_NORM_batch_UMAP.h5ad"
    
    # 1. Quality control
    print(f"Reading data: {INPUT_PATH}")
    adata = sc.read_h5ad(INPUT_PATH)
    adata = quality_control(adata)
    adata.write_h5ad(QC_PATH)
    print(f"Post-QC data: {adata.n_obs} cells, {adata.n_vars} genes.")
    
    # 2. Normalization and feature selection
    adata = normalize_and_select_features(adata)
    adata.write_h5ad(NORMALIZED_PATH)
    print(f"Normalized data saved to: {NORMALIZED_PATH}")
    
    # 3. Data scaling and PCA
    adata, n_pcs = scale_and_pca(adata)
    
    # 4. Batch correction (if batch information is available)
    if 'batch' in adata.obs.columns:
        adata = batch_correction_with_scvi(adata, batch_key='batch')
    else:
        print("Warning: No batch information found in the data, skipping batch correction")
        adata.obsm['X_scVI'] = adata.obsm['X_pca'][:, :30]  # Use PCA as a fallback
    
    adata.write_h5ad(BATCH_CORRECTED_PATH)
    print(f"Batch-corrected data saved to: {BATCH_CORRECTED_PATH}")
    
    # 5. Cell clustering and UMAP visualization
    adata = cluster_and_visualize(adata)
    adata.write_h5ad(CLUSTERED_PATH)
    print(f"Clustering and visualization results saved to: {CLUSTERED_PATH}")
    
    # 6. Generate UMAP visualization plot
    sc.pl.umap(adata, color=['louvain', 'batch'], save='_clusters.png')
    print("UMAP visualization saved.")

if __name__ == "__main__":
    main()