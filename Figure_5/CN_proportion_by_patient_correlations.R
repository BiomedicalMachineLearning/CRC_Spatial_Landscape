library(corrplot)
library(ggbeeswarm)
library(Hmisc)
library(tidyverse)

sessionInfo()

set.seed(42)

# Load data ---------------------------------------------------------------

DATA_DIR <- './CRC_Spatial_Landscape/Data/'
OUTPUT_DIR <- './CRC_Spatial_Landscape/Figure_5'

samples <- read_csv(file.path(DATA_DIR, 'samples.csv'))

cn_props <- read_csv(file.path(DATA_DIR, "cell_neighbourhoods_k8_proportions_by_patient.csv"))

# Preprocess cell neighbourhood proportions -------------------------------

## Filter for cancer samples -----------------------------------------------

cn_props <- cn_props %>%
  filter(batch=='Cancer') %>%
  left_join(samples %>% select(alt_identifier,
                               MSI_any,
                               Reccurrance),
            by='alt_identifier')
  
# Correlation functions ---------------------------------------------------

calc_cn_corr <- function(x){
  x_m <- x %>%
    select(`Bulk Tumor`:Stromal) %>%
    as.matrix()
  res <- rcorr(x_m, type='spearman')
  
  # FDR correction for multiple testing
  res$P_adj <-  res$P %>%
    p.adjust(method='fdr') %>%
    matrix(ncol = ncol(res$r))
  dimnames(res$P_adj) <- dimnames(res$r)
  return(res)
}

plot_corr <- function(res){
  corrplot(res$r, type = 'upper', order = 'alphabet', hclust.method = 'ward.D2', addrect=3,
           diag = FALSE, p.mat = res$P_adj, sig.level = 0.05, insig = 'blank', pch.cex = 0.9, tl.col='black', tl.srt = 45,
           col = palette_corr)
}

palette_corr <- rev(COL2('RdBu', 200))

# All patients ------------------------------------------------------------

pdf(file.path(OUTPUT_DIR, "CN_proportion_by_patient_correlation_fdr_all_cancers.pdf"), width = 6, height=6)
cn_props %>%
  calc_cn_corr() %>%
  plot_corr()
dev.off()

# MSI status --------------------------------------------------------------

pdf(file.path(OUTPUT_DIR, "CN_proportion_by_patient_correlation_fdr_MSI.pdf"), width = 6, height=6)
cn_props %>%
  filter(MSI_any == 'TRUE') %>%
  calc_cn_corr() %>%
  plot_corr()
dev.off()

pdf(file.path(OUTPUT_DIR, "CN_proportion_by_patient_correlation_fdr_MSS.pdf"), width = 6, height=6)
cn_props %>%
  filter(MSI_any == 'FALSE') %>%
  calc_cn_corr() %>%
  plot_corr()
dev.off()

# Recur status ------------------------------------------------------------

pdf(file.path(OUTPUT_DIR, "CN_proportion_by_patient_correlation_fdr_recur.pdf"), width = 6, height=6)
cn_props %>%
  filter(Reccurrance == 'recur') %>%
  calc_cn_corr() %>%
  plot_corr()
dev.off()

pdf(file.path(OUTPUT_DIR, "CN_proportion_by_patient_correlation_fdr_no_recur.pdf"), width = 6, height=6)
cn_props %>%
  filter(Reccurrance == 'no recur') %>%
  calc_cn_corr() %>%
  plot_corr()
dev.off()

# MSI-H Recur and No Recur -------------------------------------------------

pdf(file.path(OUTPUT_DIR, "CN_proportion_by_patient_correlation_fdr_MSI_recur.pdf"), width = 6, height=6)
cn_props %>%
  filter(MSI_any == TRUE, 
         Reccurrance == 'recur') %>%
  calc_cn_corr() %>%
  plot_corr()
dev.off()

# Not enough samples for significance testing in MSI no recur patients

# MSS Recur and No Recur ---------------------------------------------------

pdf(file.path(OUTPUT_DIR, "CN_proportion_by_patient_correlation_fdr_MSS_recur.pdf"), width = 6, height=6)
cn_props %>%
  filter(MSI_any == FALSE, 
         Reccurrance == 'recur') %>%
  calc_cn_corr() %>%
  plot_corr()
dev.off()

pdf(file.path(OUTPUT_DIR, "CN_proportion_by_patient_correlation_fdr_MSS_no_recur.pdf"), width = 6, height=6)
cn_props %>%
  filter(MSI_any == FALSE, 
         Reccurrance == 'no recur') %>%
  calc_cn_corr() %>%
  plot_corr()
dev.off()
