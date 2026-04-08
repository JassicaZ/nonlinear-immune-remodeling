args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]
date <- args[2]

# Load required libraries
library('DEswan')
library(edgeR)
library(MASS)


# Read input files
## Read clinical information
cli <- read.csv('./data/all_sample_cli_info.csv', row.names = 1)

## Read gene expression data
df <- read.csv(paste0("../result/1_preprocess/discovery_cohort_adjusted/", ct, ".csv"), row.names = 1)

## Prepare the matrix
mat <- as.matrix(df)

# Prepare data for DE-SWAN
mat4deswan <- t(mat)  # Transpose the matrix
cli4deswan <- cli[rownames(mat4deswan), ]

# Perform DE-SWAN analysis (window size = 8, sliding every 2 years)
start_time <- Sys.time()
res.DEswan <- DEswan(
  data.df = mat4deswan,
  qt = cli4deswan$age,
  window.center = seq(2, 90, 2),
  buckets.size = 8
)
end_time <- Sys.time()
end_time - start_time

# Save the results
file_path <- '../result/4_DESWAN_for_gene_expression/res/'
if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
}
saveRDS(res.DEswan, paste0(file_path, ct, '.rds'))

# Process the results
## Convert long format to wide format
res.DEswan.wide.p <- reshape.DEswan(res.DEswan, parameter = 1, factor = "qt")

## Adjust p-values
res.DEswan.wide.q <- q.DEswan(res.DEswan.wide.p, method = "BH")

# Save the table of significant variables for different thresholds
# Column names represent significance thresholds, and values represent the number of significant genes at each center point
file_path <- '../result/4_DESWAN_for_gene_expression/sig_df/'
if (!dir.exists(file_path)) {
  dir.create(file_path, recursive = TRUE)
}
write.csv(res.DEswan.wide.q.signif, paste0(file_path, ct, '.csv'))


