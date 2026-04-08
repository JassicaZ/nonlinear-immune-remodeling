library(rlang)
library(Mfuzz)
library(dplyr)
library(stringr)
library(data.table)


# Get cell type
ctfilepath <- '../result/1_preprocess/discovery_cohort_adjusted/sig_df/'
files <- list.files(path = ctfilepath, pattern = '*.csv', full.names = TRUE)
cts <- sub(paste0(".*/(.*)\\.csv$"), "\\1", files)

# Merge all tables with gene expression tradjories
all_df <- data.frame()
for (ct in cts) {
    if (file.exists(paste0("../result/3_trajectory/loess/", ct, ".csv"))) {
        df <- data.frame(fread(paste0("../result/3_trajectory/loess/", ct, ".csv")))
        colnames(df)[1] <- 'gene'
        df$gene <- paste0(ct, '_', df$gene)
        all_df <- rbind(all_df, df[-1, ])
    }
}

rownames(all_df) <- all_df$gene
all_df <- all_df[-1]
colnames(all_df) <- 1:100

# Merge all significant nonlinear genes
sig_df <- data.frame()
for (ct in cts) {
    df <- data.frame(fread(paste0("../result/2_age_associated_genes/non_linear_association/sig_df/", ct, ".csv")))
    sig_gene <- rownames(df)
    colnames(df)[1] <- 'gene'
    df$gene <- paste0(ct, '_', df$gene)
    sig_df <- rbind(sig_df, df[-1, ])
}

# Filter the genes significant with age
selected_df <- all_df[rownames(all_df) %in% sig_df$gene, ]

# Perform Mfuzz clustering analysis
mat <- as.matrix(t(scale(t(selected_df))))
eset <- Biobase::ExpressionSet(mat)

eset <- standardise(eset)

m <- mestimate(eset)

dim <- Mfuzz::Dmin(eset, m, crange = seq(2, 20, 2), repeats = 1, visu = TRUE)

result <- mfuzz(eset, centers = 6, m = m)

# Generate clustering plots
plot_mfuzz <- function(eset, result, group_name) {
    pdf(paste0("../result/3_trajectory/mfuzz/", group_name, '_mfuzz_clusters.pdf'), width = 10, height = 12)
    mfuzz.plot2(eset, cl = result, x11 = FALSE, mfrow = c(3, 2), centre = TRUE, xlab = "Age (years)", ylab = "z-score")
    dev.off()
}

# Generate clustering plots
plot_mfuzz(eset, result, 'sig_6cl')