library('DEswan')
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)


# Process metabolomics data
# Retain metabolites present in the HMDB or KEGG databases
info_df <- read.csv('../data/bulk_serum/metabolite_info.csv')
info_df_filtered <- info_df[!is.na(info_df$HMDB.ID) | !is.na(info_df$KEGG.ID), ]

# Prepare metabolomics data matrix and corresponding clinical information
metabo_df <- read.csv('../data/bulk_serum/metabonomics_serum.csv', row.names = 1)
cli <- read.csv('../data/bulk_serum/cli_age_sex.csv', row.names = 1)
ID_df <- read.csv('../data/bulk_serum/ID_serum_file.csv')

metabo_df_filtered <- metabo_df[
  rownames(metabo_df) %in% info_df_filtered$Name,
  colnames(metabo_df) %in% ID_df$GBI_serum_id
]

rownames(ID_df) <- ID_df[, 1]
ID_df <- ID_df[, -1, drop = FALSE]

metabo_df_filtered <- metabo_df_filtered[, cli$GBI_serum_id]

# Run DE-SWAN
mat <- as.matrix(metabo_df_filtered)
cli <- na.omit(cli)
rownames(cli) <- cli$GBI_serum_id

mat4deswan <- t(log10(mat + 1))
cli4deswan <- cli[rownames(mat4deswan), ]

start_time <- Sys.time()
res.DEswan <- DEswan(
  data.df = mat4deswan,
  qt = cli4deswan$age,
  window.center = seq(4, 20, 2),
  buckets.size = 4,
  covariates = cli4deswan['sex']
)
end_time <- Sys.time()
end_time - start_time

# Process results
res.DEswan.wide.p <- reshape.DEswan(res.DEswan, parameter = 1, factor = "qt")
res.DEswan.wide.q <- q.DEswan(res.DEswan.wide.p, method = "BH")
res.DEswan.wide.q.signif <- nsignif.DEswan(res.DEswan.wide.q)

saveRDS(res.DEswan, '../result/4_DESWAN_for_metabolites_and_lipids/res/metabolites_win4_len2_res.rds')
write.csv(res.DEswan.wide.q.signif, '../result/4_DESWAN_for_metabolites_and_lipids/sig_df/metabolites_win4_len2_sig.csv')

# Process lipidomics data
lipid_df <- read.csv('../data/bulk_serum/data_Lipid_HC.csv', row.names = 1)
cli <- read.csv('../data/bulk_serum/cli_age_sex.csv', row.names = 1)
lipid_df <- lipid_df[, cli$GBI_serum_id]

start_time <- Sys.time()
res.DEswan <- DEswan(
  data.df = mat4deswan,
  qt = cli4deswan$age,
  window.center = seq(4, 20, 2),
  buckets.size = 4,
  covariates = cli4deswan['sex']
)
end_time <- Sys.time()
end_time - start_time

# Process results
res.DEswan.wide.p <- reshape.DEswan(res.DEswan, parameter = 1, factor = "qt")
res.DEswan.wide.q <- q.DEswan(res.DEswan.wide.p, method = "BH")
res.DEswan.wide.q.signif <- nsignif.DEswan(res.DEswan.wide.q)


saveRDS(res.DEswan, '../result/4_DESWAN_for_metabolites_and_lipids/res/lipids_win4_len2_res.rds')
write.csv(res.DEswan.wide.q.signif, '../result/4_DESWAN_for_metabolites_and_lipids/sig_df/lipids_win4_len2_sig.csv')