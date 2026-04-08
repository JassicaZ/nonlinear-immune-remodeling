library(Hmisc)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)

#Association analysis between gene expression and metabolite levels
sig_gene_df <- read.csv('../result/4_DESWAN_for_gene_expression/sig_genes_in_each_age.csv',row.names=1)
sig_gene_df <- sig_gene_df[sig_gene_df$age %in% c('X6','X8'), ]
sig_meta_df <- read.csv('../result/5_DESWAN_for_metabolites_and_lipids/sig_df/metabolite_age68_sig_with_info.csv',row.names=1)

meta_df <- read.csv('../data/bulk_serum/metabonomics_serum.csv',row.names=1)
meta_df <- meta_df[sig_meta_df$variable,cli_df$GBI_serum_id]
colnames(meta_df) <- cli_df$RNA_id
M <- log10(meta_df+1)

genes <- rownames(all_df)
chunk_size <- 500
df_p <-data.frame()
df_r <-data.frame()
n_meta <- dim(M)[1]
for (i in seq(1, length(genes), by = chunk_size)) {
  chunk_genes <- genes[i:min(i + chunk_size - 1, length(genes))]
  sub_mat <- all_df[chunk_genes, ]
  sub_res <- rcorr(t(M), t(sub_mat), type="spearman")
  
  #Row for genes, column for metabolites
  row_idx <- (n_meta + 1):(n_meta + length(chunk_genes))
  col_idx <- 1:n_meta
  
  df_p <- rbind(df_p, sub_res$P[row_idx, col_idx])
  df_r <- rbind(df_r, sub_res$r[row_idx, col_idx]) 

}

df_p <- df_p[,!grepl('\\.1',colnames(df_p))]
df_r <- df_r[,!grepl('\\.1',colnames(df_r))]

df_p_long <- df_p %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(
    cols = colnames(df_p),
    names_to = "metabolite",
    values_to = "p_value"
  )%>%
  group_by(gene) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

  df_r_long <- df_r %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(
    cols = colnames(df_r),
    names_to = "metabolite",
    values_to = "cor"
  )


df_p_r <-merge(df_p_long,df_r_long,by=c('gene','metabolite'))
df_p_r<- df_p_r[df_p_r$p_adj<0.05 & abs(df_p_r$cor)>0.2,]

df_p_r$pair <- paste0(df_p_r$celltype,'_',df_p_r$gene)
df_p_r_X68<- df_p_r[df_p_r$pair %in% sig_gene_df$pair,]

write.csv(df_p_r_X68,'../result/5_DESWAN_for_metabolites_and_lipids/rna_meta_spearman_age68.csv')


#Association analysis between gene expression and lipid levels

sig_gene_df <- read.csv('../result/4_DESWAN_for_gene_expression/sig_genes_in_each_age.csv',row.names=1)
sig_gene_df <- sig_gene_df[sig_gene_df$age == 'X6', ]
sig_meta_df <- read.csv('../result/5_DESWAN_for_metabolites_and_lipids/sig_df/lipid_age6_sig_with_info.csv',row.names=1)

meta_df <- read.csv('../data/bulk_serum/metabonomics_serum.csv',row.names=1)
meta_df <- meta_df[sig_meta_df$variable,cli_df$GBI_serum_id]
colnames(meta_df) <- cli_df$RNA_id
M <- log10(meta_df+1)

genes <- rownames(all_df)
chunk_size <- 500
df_p <-data.frame()
df_r <-data.frame()
n_meta <- dim(M)[1]
for (i in seq(1, length(genes), by = chunk_size)) {
  chunk_genes <- genes[i:min(i + chunk_size - 1, length(genes))]
  sub_mat <- all_df[chunk_genes, ]
  sub_res <- rcorr(t(M), t(sub_mat), type="spearman")
  
  #Row for genes, column for metabolites
  row_idx <- (n_meta + 1):(n_meta + length(chunk_genes))
  col_idx <- 1:n_meta
  
  df_p <- rbind(df_p, sub_res$P[row_idx, col_idx])
  df_r <- rbind(df_r, sub_res$r[row_idx, col_idx]) 

}

df_p <- df_p[,!grepl('\\.1',colnames(df_p))]
df_r <- df_r[,!grepl('\\.1',colnames(df_r))]

df_p_long <- df_p %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(
    cols = colnames(df_p),
    names_to = "lipid",
    values_to = "p_value"
  )%>%
  group_by(gene) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

  df_r_long <- df_r %>%
  rownames_to_column(var = "gene") %>% 
  pivot_longer(
    cols = colnames(df_r),
    names_to = "lipid",
    values_to = "cor"
  )


df_p_r <-merge(df_p_long,df_r_long,by=c('gene','lipid'))
df_p_r<- df_p_r[df_p_r$p_adj<0.05 & abs(df_p_r$cor)>0.2,]

df_p_r$pair <- paste0(df_p_r$celltype,'_',df_p_r$gene)
df_p_r_X6<- df_p_r[df_p_r$pair %in% sig_gene_df$pair,]

write.csv(df_p_r_X6,'../result/5_DESWAN_for_metabolites_and_lipids/rna_lipid_spearman_age6.csv')