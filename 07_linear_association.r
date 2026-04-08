library('Hmisc')
library(dplyr)
library(stringr)

# Get cell types
ctfilepath<- '../data/1_preprocess/discovery_cohort_adjusted/'
files <- list.files(path=ctfilepath, pattern='*.csv',full.names=TRUE)
cts <- sub(paste0(".*/(.*)\\.csv$"), "\\1", files)

for (ct in cts) {
print(ct)

    cli<- read.csv('../data/discovery_sample_cli_info.csv', row.names = 1)

    residuals_df <- read.csv(paste0("../result/1_preprocess/discovery_cohort_adjusted/",ct,".csv"), row.names = 1)
    residuals_mat <- as.matrix(residuals_df)

    cli<-cli[colnames(residuals_mat),]
    age_vector <- as.vector(cli$age)
#Perform linear correlation analysis
results <- apply(residuals_mat, 1, function(x) {
  ct <- cor.test(age_vector, x, method = "pearson")
  c(cor = ct$estimate, p = ct$p.value)
})
df <- data.frame(t(results))   
df$p.adj <- p.adjust(df$p, method='BH')
write.csv(df,paste0('./result/2_age_associated_genes/linear_association/',ct,'.csv'))

df005 <- df[abs(df$cor.cor)>0.2 & df$p.adj <0.05,]
write.csv(df005,paste0('./result/2_age_associated_genes/linear_association/padj005/',ct,'.csv'))

}