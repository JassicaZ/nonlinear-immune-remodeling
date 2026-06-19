import scanpy as sc
import numpy as np
import os
import celltypist
from celltypist import models
import time  

# Set the folder path
file_path = '../result/1_preprocess/celltypist/'

date = '250520_customized'
folder_path = file_path + '/{}'.format(date)
if not os.path.exists(folder_path):
    # Create the folder if it does not exist
    os.makedirs(folder_path)
    print("Folder has been created.")
else:
    print("Folder already exists.")
sc.settings.figdir = "../result/1_preprocess/{}/".format(date)

# Model training
adata_aida = sc.read('./data/AIDA_raw_with-meta.h5ad')  # From Kock, Kian Hong et al. Cell 2025, available at https://cellxgene.cziscience.com/collections/ced320a1-29f3-47c1-a735-513c7084d508
adata_aida.var.index = adata_aida.var['feature_name'].str.replace(r'_ENSG.*$', '', regex=True)
sc.pp.normalize_total(adata_aida, target_sum=10**4)  
sc.pp.log1p(adata_aida)
adata_aida.var.index.name = None

# Downsample the data
sampled_cell_index = celltypist.samples.downsample_adata(
    adata_aida, mode='each', n_cells=1000, by='Annotation_Level3', return_index=True
)
model_fs = celltypist.train(
    adata_aida[sampled_cell_index], 'Annotation_Level3', n_jobs=10, max_iter=5, use_SGD=True
)

gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -300, axis = 1)[:, -300:]
gene_index = np.unique(gene_index)

adata4model = adata_aida[sampled_cell_index, gene_index]

# Train the model
model = celltypist.train(
    adata4model, 'Annotation_Level3', check_expression=False, n_jobs=10, max_iter=200
)

# Save the model
print(f"Time elapsed: {(t_end - t_start)/60} minutes")
model.write('./data/train/model_from_AIDA_03_v2.pkl')

# Load the trained model
model_3 = models.Model.load(model="./data/train/model_from_AIDA_03_v2.pkl")

# Preprocess discovery cohort single-cell data
adata_raw = sc.read('../result/1_preprocess/preprocessed/inner_QC.h5ad')
adata_celltypist = adata_raw.copy()  # Make a copy of the AnnData object
adata_celltypist.X = adata_raw.layers["counts"]  # Set adata.X to raw counts
sc.pp.normalize_total(adata_celltypist, target_sum=10**4)  # Normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # Log-transform the data
adata_celltypist.X = adata_celltypist.X.toarray()  # Convert sparse matrix to dense for compatibility with celltypist
adata_celltypist.write(file_path + '/{}/{}_adata_celltypist.h5ad'.format(date, date))

# Load neighbors
adata = sc.read('../result/1_preprocess/preprocessed/harmony_inner.h5ad')
adata_celltypist.uns['neighbors'] = adata.uns['neighbors']
adata_celltypist.obsp['connectivities'] = adata.obsp['connectivities']
adata_celltypist.obsp['distances'] = adata.obsp['distances']

# Predict cell types
print("Performing annotation using the trained model.")
start_time = time.time()  # Start timing
print("Start time: ", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time)))  # Print start time
predictions = celltypist.annotate(
    adata_celltypist, model=model_3, majority_voting=False, use_GPU=False
)
predictions_adata = predictions.to_adata()

# Save the annotated data
predictions_adata.write(file_path + '/{}/{}_celltypist_annotated_3.h5ad'.format(date, date))

end_time = time.time()  # End timing
print("End time: ", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time)))  # Print end time
print(f"Annotation computation time: {end_time - start_time} seconds")  # Print elapsed time

# Preprocess validation cohort single-cell data
adata_raw = sc.read('../result/1_preprocess/validation/PJ_QC.h5ad')
adata_celltypist = adata_raw.copy()  # Make a copy of the AnnData object
del adata_raw
sc.pp.normalize_total(adata_celltypist, target_sum=10**4)  # Normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # Log-transform the data
adata_celltypist.X = adata_celltypist.X.toarray()  # Convert sparse matrix to dense for compatibility with celltypist
adata_celltypist.write(file_path + '/{}/{}_adata_celltypist.h5ad'.format(date, date))

# Load neighbors from another AnnData object
adata = sc.read('../result/1_preprocess/validation/PJ_QC_NORM_batch_UMAP.h5ad')
adata_celltypist.uns['neighbors'] = adata.uns['neighbors']
adata_celltypist.obsp['connectivities'] = adata.obsp['connectivities']
adata_celltypist.obsp['distances'] = adata.obsp['distances']

# Predict cell types
predictions = celltypist.annotate(
    adata_celltypist, model=model_3, majority_voting=False, use_GPU=False
)
predictions_adata = predictions.to_adata()

# Save the annotated data
predictions_adata.write(file_path + '/{}/{}_celltypist_annotated_3_30_1.3_pu.h5ad'.format(date, date))
