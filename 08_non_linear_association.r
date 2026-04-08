library(dplyr)
library(stringr)
library(data.table)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(dplyr)
library(stringr)
library(data.table)

# Limit the number of threads for a single task to prevent overloading
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]

# Read data
df <- read.csv(paste0("../result/1_preprocess/discovery_cohort_adjusted/", ct, ".csv"), row.names = 1)
mat <- as.matrix(df)

counts <- mat
cli <- read.csv('../data/discovery_sample_cli_info.csv', row.names = 1)
cli <- cli[colnames(counts), ]  # Only include data for individuals present in the dataset
pseudotime <- as.matrix(cli$age)  # Use age as pseudotime

# Construct cellWeights, all set to 1 (since each cell follows only one "age trajectory")
cellWeights <- matrix(1, nrow = nrow(pseudotime), ncol = 1)
colnames(cellWeights) <- "age_lineage"
colnames(pseudotime) <- "age_lineage"

# Covariates
cli$pub <- factor(cli$pub)
cli$sex <- factor(cli$sex)

offset <- rep(log(1), times = dim(mat)[2])

# Set up parallel computing
# Modify fitGAM source code for parallel computing, see https://github.com/statOmics/tradeSeq/issues/261
source('fitGAM_modified.R')
source('_fitGAM_modified.R')

library(batchtools)
library(BiocParallel)
library(tictoc) 

param <- BatchtoolsParam()

tic("fitGAM modified")

# Get the total number of genes
total_genes <- nrow(counts)
cat("Total genes:", total_genes, "\n")

# Set batch size (adjust based on CPU limits, e.g., process 50-100 genes per batch)
batch_size <- 30  # Adjust based on your CPU limits

# Process in batches
results_list <- list()

for (i in seq(1, total_genes, batch_size)) {
  end_idx <- min(i + batch_size - 1, total_genes)
  cat(sprintf("Processing batch: genes %d to %d (%d/%d)\n", 
              i, end_idx, ceiling(i / batch_size), ceiling(total_genes / batch_size)))
  print(end_idx)
  
  # Extract the current batch of genes
  batch_expr <- counts[i:end_idx, , drop = FALSE]
  batch_expr <- t(scale(t(batch_expr), center = TRUE, scale = TRUE))  # Scale gene expression to prevent prediction overflow
  
  # Run fitGAM on the current batch
  batch_result <- fitGAM_modified(
    counts = batch_expr,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    nknots = 5,
    offset = offset,
    parallel = TRUE,
    BPPARAM = param,
    verbose = TRUE,
    family = "gaussian"
  )
  
  # Save the results of the current batch
  results_list[[length(results_list) + 1]] <- batch_result
  
  cat(sprintf("Batch %d completed\n", ceiling(i / batch_size)))
}

# Combine results from all batches
cat("Combining results from all batches...\n")
sce <- do.call(rbind, results_list)

toc()


# Perform association testing
rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, contrastType = "consecutive")
rowData(sce)$assocRes$p.adj <- p.adjust(rowData(sce)$assocRes$pvalue_1, method = 'BH')

# Save the results
saveRDS(sce, paste0("../result/2_age_associated_genes/non_linear_association/sce/", ct, ".rds"))

assoRes <- rowData(sce)$assocRes
sig_df <- assoRes[!is.na(assoRes$p.adj) & assoRes$p.adj < 0.05, ]
sig_gene <- rownames(sig_df)

filename <- paste0("../result/2_age_associated_genes/non_linear_association/sig_df/", ct, ".csv")
write.csv(sig_df, file = filename)

