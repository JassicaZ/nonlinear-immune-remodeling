library(gamlss)
library(gamlss.add)
library(gamlss.dist)
library(purrr)
library(tidyverse)
library(ggbeeswarm)
library(gghalves)
library(glue)

input_dir <- '../result/5_peak_gene_model/'

# 1. Load training and validation data
train_df <- read.csv(paste0(input_dir, 'train_result.csv'), row.names = 1)

cohort_names <- c('SLE', 'CD', 'T1D', 'test')
cohort_list <- list()

for (cohort in cohort_names) {
  cohort_list[[cohort]] <- read.csv(paste0(input_dir, cohort, '_result.csv'))
}

# 2. Model selection: compare LMS, LMSP, and LMST families
## 2.1 Fit models
fit_LMS <- gamlss(
  remodeling_score ~ pb(age),
  sigma.formula = ~ pb(age),
  nu.formula = ~ pb(age),
  family = BCCG,
  data = train_df,
  control = gamlss.control(n.cyc = 100, trace = FALSE)
)

fit_LMSP <- gamlss(
  remodeling_score ~ pb(age),
  sigma.formula = ~ pb(age),
  nu.formula = ~ pb(age),
  tau.formula = ~ pb(age),
  family = BCPE,
  data = train_df,
  control = gamlss.control(n.cyc = 100, trace = FALSE)
)

fit_LMST <- gamlss(
  remodeling_score ~ pb(age),
  sigma.formula = ~ pb(age),
  nu.formula = ~ pb(age),
  tau.formula = ~ pb(age),
  family = BCT,
  data = train_df,
  control = gamlss.control(n.cyc = 100, trace = FALSE)
)

## 2.2 Compare models using AIC and BIC
AIC(fit_LMS, fit_LMSP, fit_LMST)
BIC(fit_LMS, fit_LMSP, fit_LMST)

## 2.3 Visual diagnostics for each model
centiles_values <- c(5, 10, 25, 50, 75, 90, 95)

centiles(fit_LMS, xvar = train_df$age, cent = centiles_values, main = 'LMS (BCCG)')
centiles(fit_LMSP, xvar = train_df$age, cent = centiles_values, main = 'LMSP (BCPE)')
centiles(fit_LMST, xvar = train_df$age, cent = centiles_values, main = 'LMST (BCT)')

options(repr.plot.width = 8, repr.plot.height = 6)

wp(fit_LMS, xvar = train_df$age, n.inter = 9, main = 'Worm plot - LMS')
wp(fit_LMSP, xvar = train_df$age, n.inter = 9, main = 'Worm plot - LMSP')
wp(fit_LMST, xvar = train_df$age, n.inter = 9, main = 'Worm plot - LMST')

plot(fit_LMS)
plot(fit_LMSP)
plot(fit_LMST)
#Final model selection: BCPE (LMSP)

# 3. model fitting with BCPE family (LMSP)
fit_LMSP <- gamlss(
  remodeling_score ~ pb(age),
  sigma.formula = ~ pb(age),
  nu.formula = ~ pb(age),
  tau.formula = ~ pb(age),
  family = BCPE,
  data = train_df,
  control = gamlss.control(n.cyc = 100, trace = FALSE)
)
saveRDS(fit_LMSP, paste0(input_dir, 'fit_LMSP.rds'))


# 4. Z-score calculation
## 4.1 Function to calculate Z-score based on BCPE model
get_zscore_BCPE <- function(model, newdata) {
  newdata_pred <- newdata[, 'age', drop = FALSE]
  pred <- predictAll(model, newdata = newdata_pred, type = 'response', data = train_df)

  y <- newdata$remodeling_score
  mu <- pred$mu
  sigma <- pred$sigma
  nu <- pred$nu
  tau <- pred$tau

  t <- ifelse(
    abs(nu) > 1e-6,
    ((y / mu)^nu - 1) / (nu * sigma),
    log(y / mu) / sigma
  )

  p <- pPE(t, mu = 0, sigma = 1, nu = tau)
  qnorm(p)
}


## 4.2 Z-score for training cohort 
train_zscores <- get_zscore_BCPE(fit_LMSP, train_df)

train_summary <- tibble(
  diagnosis = 'train',
  age = train_df$age,
  remodeling_score = train_df$remodeling_score,
  zscore = train_zscores,
  display_label = 'Train'
)

## 4.3 Z-score for validation cohorts
DIAGNOSIS_LABELS <- c(
  SLE = 'SLE',
  CD = 'CD',
  T1D = 'T1D',
  test = 'Test'
)

zscore_long <- imap_dfr(cohort_list, function(df, disease_name) {
  message('calculate z-score: ', disease_name)

  df_clean <- df %>% filter(!is.na(age), !is.na(remodeling_score))
  z <- get_zscore_BCPE(fit_LMSP, df_clean)

  tibble(
    diagnosis = disease_name,
    age = df_clean$age,
    remodeling_score = df_clean$remodeling_score,
    zscore = z
  )
}) %>%
  mutate(
    display_label = factor(
      DIAGNOSIS_LABELS[diagnosis],
      levels = DIAGNOSIS_LABELS[names(cohort_list)]
    )
  )

zscore_long %>%
  group_by(diagnosis) %>%
  summarise(
    n = n(),
    mean = round(mean(zscore, na.rm = TRUE), 3),
    sd = round(sd(zscore, na.rm = TRUE), 3),
    median = round(median(zscore, na.rm = TRUE), 3),
    q25 = round(quantile(zscore, 0.25, na.rm = TRUE), 3),
    q75 = round(quantile(zscore, 0.75, na.rm = TRUE), 3),
    .groups = 'drop'
  ) %>%
  print()

# 5. Permutation test for significance
N_PERM <- 10000
SET_SEED <- 42
SIG_ALPHA <- 0.05

train_zscores <- train_summary %>% pull(zscore)
n_train <- length(train_zscores)

run_permutation <- function(disease_zscore, ref_zscore = train_zscores, n_perm = N_PERM) {
  n_disease <- length(disease_zscore)
  observed_diff <- median(disease_zscore) - median(ref_zscore)
  combined <- c(ref_zscore, disease_zscore)
  n_total <- length(combined)

  perm_diffs <- replicate(n_perm, {
    shuffled <- sample(combined)
    median(shuffled[(n_train + 1):n_total]) - median(shuffled[1:n_train])
  })

  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))

  list(
    n_disease = n_disease,
    observed_diff = observed_diff,
    p_value = p_value,
    perm_diffs = perm_diffs
  )
}

perm_results <- map_dfr(names(cohort_list), function(dname) {
  message('permutation test: ', dname)

  disease_zscore <- zscore_long %>%
    filter(diagnosis == dname) %>%
    pull(zscore)

  res <- run_permutation(disease_zscore)

  tibble(
    diagnosis = dname,
    n = length(disease_zscore),
    median_disease = median(disease_zscore),
    median_test = median(train_zscores),
    observed_diff = res$observed_diff,
    p_value = res$p_value
  )
}) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = 'BH'),
    significance = case_when(
      p_adjusted < 0.001 ~ '***',
      p_adjusted < 0.01 ~ '**',
      p_adjusted < 0.05 ~ '*',
      TRUE ~ 'ns'
    )
  )

message('\nPermutation test results (BH-FDR adjusted):')

perm_results %>%
  select(diagnosis, n, median_disease, observed_diff, p_value, p_adjusted, significance) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print()
