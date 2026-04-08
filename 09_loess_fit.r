library(stats)
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]

# Define a function to select the optimal span for LOESS using cross-validation
cv_loess_fast <- function(x, y, spans = seq(0.2, 0.6, 0.1), k = 5) {
  n <- length(y)
  folds <- sample(rep(1:k, length.out = n))  # Create k-folds for cross-validation
  
  # Compute RMSE for each span
  cv_rmse <- sapply(spans, function(span) {
    errs <- numeric(k)
    for (i in 1:k) {
      test_idx <- which(folds == i)
      fit <- loess(y ~ x, data = data.frame(x, y), subset = setdiff(1:n, test_idx), span = span)
      pred <- predict(fit, newdata = data.frame(x = x[test_idx]))
      errs[i] <- mean((y[test_idx] - pred)^2, na.rm = TRUE)
    }
    sqrt(mean(errs))  # Return RMSE for the current span
  })
  
  # Select the span with the minimum RMSE
  best_span <- spans[which.min(cv_rmse)]
  return(best_span)
}

# Read input data
df <- read.csv(paste0("../result/1_preprocess/discovery_cohort_adjusted/", ct, ".csv"), row.names = 1)
mat <- as.matrix(df)

cli <- read.csv('../data/discovery_sample_cli_info.csv', row.names = 1)
cli <- cli[colnames(mat), ]

# Initialize variables
gene_n <- dim(mat)[1]
all_df <- data.frame()

# Perform LOESS fitting for each gene
for (i in 1:gene_n) {
  y <- mat[i, ]
  x <- cli$age
  
  span <- cv_loess_fast(x, y, spans = seq(0.2, 1, 0.1), k = 5)
  
  # Fit LOESS model and predict values
  fit <- loess(y ~ x, span = span)
  pred <- predict(fit, newdata = data.frame(x = seq(min(x), max(x), length.out = 100)))
  
  pred_df <- t(data.frame(pred))
  rownames(pred_df) <- rownames(mat)[i]
  all_df <- rbind(all_df, pred_df)
  
  print(paste0('Calculating ', i, ' of ', gene_n))
}

# Save the results
output_path <- paste0("../result/3_trajectory/loess/", ct, ".csv")
write.csv(all_df, output_path)

print('LOESS fitting completed.')