library(lme4)
library(dplyr)
library(stringr)
library(data.table)
library(parallel)
library(edgeR)

## Discovery cohort

# Get cell types
ctfilepath <- '../result/1_preprocess/pseudobulk/discovery_cohort/'
files <- list.files(path = ctfilepath, pattern = '*.csv', full.names = TRUE)
cts <- sub(paste0(".*/(.*)\\.csv$"), "\\1", files)

for (ct in cts) {
    # Read data
    cli_df <- read.csv('../data/all_sample_cli_info_filtr.csv', row.names = 1)
    df <- read.csv(paste0("../result/1_preprocess/pseudobulk/discovery_cohort/", ct, ".csv"), row.names = 1)
    mat <- as.matrix(df)
    cli <- cli_df[rownames(df), ]
    age_vector <- as.vector(cli$age)
    
    # Data cleaning and TMM normalization
    mat <- mat[, colSums(mat != 0) > 100]  # Keep genes expressed in at least 100 individuals
    group <- rownames(mat)
    d <- DGEList(counts = t(mat), group = group)
    TMM <- calcNormFactors(d, method = "TMM")
    mat <- cpm(TMM, log = TRUE, prior.count = 1)
    
    cli$sex <- factor(cli$sex)
    cli$batch <- factor(cli$pub)
    
    # Parallel computation to obtain residuals
    batch_size <- 1000
    n_genes <- nrow(mat)
    residuals_mat <- matrix(NA, nrow = n_genes, ncol = nrow(cli))
    rownames(residuals_mat) <- rownames(mat)
    colnames(residuals_mat) <- rownames(cli)
    n_batches <- ceiling(n_genes / batch_size)
    
    for (b in 1:n_batches) {
        idx_start <- (b - 1) * batch_size + 1
        idx_end <- min(b * batch_size, n_genes)
        idx <- idx_start:idx_end
        h
        n_cores <- min(detectCores() - 1, length(idx))
        batch_res <- mclapply(idx, function(i) {
            y <- as.numeric(mat[i, ])
            cli$expr <- y
            fit <- tryCatch(lmer(expr ~ age + sex + (1 | pub), data = cli),
                            error = function(e) NULL)
            if (is.null(fit)) return(rep(NA, nrow(cli)))
            # Retain age fixed effect, remove sex and pub
            resid(fit) + fixef(fit)["(Intercept)"] + fixef(fit)["age"] * cli$age
        }, mc.cores = n_cores)
        
        # Write results to matrix
        residuals_mat[idx, ] <- do.call(rbind, batch_res)
        
        cat("Batch", b, "processed:", idx_start, "-", idx_end, "\n")
    }
    
    rownames(residuals_mat) <- rownames(mat)
    write.csv(residuals_mat, paste0('../result/1_preprocess/discovery_cohort_adjusted/', ct, '.csv'))
}

## Validation cohort
read_fread_with_rownames <- function(file, rowname_col = 1) {
  dt <- fread(file, header = TRUE)
  df <- as.data.frame(dt)
  rn <- df[[rowname_col]]
  df <- df[, -rowname_col, drop = FALSE]
  rownames(df) <- rn
  return(df)
}

for (ct in cts) { 
    cli<- read.csv('../data/all_sample_cli_info.csv', row.names = 1)
    df <- read_fread_with_rownames(paste0("../result/1_preprocess/pseudobulk/validation_cohort/",ct,".csv"))
    mat <- as.matrix(df)
    cli<-cli[rownames(df),]
    age_vector <- as.vector(cli$age)
    #Data cleaning and TMM normalization
    mat <- mat[,colSums(mat!=0)>30]#至少在20个人中有表达
    group <- rownames(mat)
    d <- DGEList(counts = t(mat), group=group)
    TMM <- calcNormFactors(d, method="TMM") 
    mat <- cpm(TMM, log = TRUE, prior.count = 1)

    cli$sex <- factor(cli$sex)
    cli$pub <- factor(rownames(cli))

# Parallel computation to obtain residuals
    batch_size <- 1000
    n_genes <- nrow(mat)
    residuals_mat <- matrix(NA, nrow = n_genes, ncol = nrow(cli))
    rownames(residuals_mat) <- rownames(mat)
    colnames(residuals_mat) <- rownames(cli)

    n_batches <- ceiling(n_genes / batch_size)

for (b in 1:n_batches) {
  idx_start <- (b - 1) * batch_size + 1
  idx_end <- min(b * batch_size, n_genes)
  idx <- idx_start:idx_end
  
  n_cores <- min(detectCores() - 1, length(idx))
  batch_res <- mclapply(idx, function(i) {
    y <- as.numeric(mat[i, ])
    cli$expr <- y
    
    fit <- tryCatch(lm(expr ~ age + sex , data = cli),
                    error = function(e) NULL)
    if (is.null(fit)) return(rep(NA, nrow(cli)))
    
    # Retain age fixed effect, remove sex
    resid(fit) + coef(fit)["(Intercept)"] + coef(fit)["age"] * cli$age
  }, mc.cores = n_cores)
  
  residuals_mat[idx, ] <- do.call(rbind, batch_res)
  
  cat("Batch", b, "processed:", idx_start, "-", idx_end, "\n")
}
rownames(residuals_mat) <- rownames(mat)
write.csv(residuals_mat,paste0('../result/1_preprocess/validation_cohort_adjusted/',ct,'.csv'))
        }

