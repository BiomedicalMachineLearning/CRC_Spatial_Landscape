library(ggbeeswarm)
library(tidyverse)

sessionInfo()

set.seed(42)

# Load data ---------------------------------------------------------------

DATA_DIR <- './CRC_Spatial_Landscape/Data/'
OUTPUT_DIR <- './CRC_Spatial_Landscape/Analysis/Cell_neighbourhoods'

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

## Transform to long format ------------------------------------------------

cn_props_long <- cn_props %>%
  pivot_longer(`Bulk Tumor`:Stromal, names_to = 'CN', values_to = 'proportion')

## Order CNs by average abundance ------------------------------------------

cn_order_sorted <- cn_props_long %>% 
  group_by(CN) %>%
  summarise(median = median(proportion)) %>%
  arrange(desc(median)) %>%
  pull(CN)

cn_props_long <- cn_props_long %>%
  mutate(CN = factor(CN, levels = cn_order_sorted))
  
# MSI-H vs MSS ------------------------------------------------------------

# Filter to only samples with MSI data
cn_props_long_msi <- cn_props_long %>%
  drop_na(MSI_any) %>%
  mutate(MSI = replace(MSI_any, MSI_any == TRUE, 'MSI-H'),
         MSI = replace(MSI, MSI == FALSE, 'MSS'),
         MSI = factor(MSI, levels = c('MSI-H', 'MSS')))

# CN proportion by MSI status
cn_props_long_msi %>%
  ggplot(mapping = aes(x=CN, y=proportion, col=MSI, fill=MSI)) +
  geom_boxplot(notch=FALSE, width=0.75, alpha=0.6, fill='white') +
  geom_quasirandom(dodge.width = 0.75, pch=21, col='black') +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
        legend.position = 'top') +
  labs(x='Cell Neighborhood', y='Proportion by Patient', col='MSI Status', fill='MSI Status')
ggsave(file.path(OUTPUT_DIR, "CN_proportion_by_patient_boxplot_MSI_vs_MSS.pdf"),
       width = 7, height=6)

# Recur vs No Recur -------------------------------------------------------

# Filter to only samples with Recur data
cn_props_long_recur <- cn_props_long %>%
  filter(Reccurrance %in% c('no recur', 'recur')) %>%
  mutate(Recur_status = replace(Reccurrance, Reccurrance == 'no recur', 'No Recur'),
         Recur_status = replace(Recur_status, Recur_status == 'recur', 'Recur'), 
         Recur_status = factor(Recur_status, levels = c('Recur', 'No Recur')))

# CN proportion by recur status
cn_props_long_recur %>%
  ggplot(mapping = aes(x=CN, y=proportion, col=Recur_status, fill=Recur_status)) +
  geom_boxplot(notch=FALSE, width=0.75, alpha=0.6, fill='white') +
  geom_quasirandom(dodge.width = 0.75, pch=21, col='black') +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
        legend.position = 'top') +
  labs(x='Cell Neighborhood', y='Proportion by Patient', col='Recur Status', fill='Recur Status')
ggsave(file.path(OUTPUT_DIR, "CN_proportion_by_patient_boxplot_recur_vs_no_recur.pdf"),
       width = 7, height=6)


# MSI-H Recur vs No Recur -------------------------------------------------

# Filter to only samples with both MSI and Recur data
cn_props_long_msi_recur <- cn_props_long_msi %>%
  filter(Reccurrance %in% c('no recur', 'recur')) %>%
  mutate(Recur_status = replace(Reccurrance, Reccurrance == 'no recur', 'No Recur'),
         Recur_status = replace(Recur_status, Recur_status == 'recur', 'Recur'), 
         Recur_status = factor(Recur_status, levels = c('Recur', 'No Recur')))

cn_props_long_msi_recur %>%
  filter(MSI == 'MSI-H') %>%
  ggplot(mapping = aes(x=CN, y=proportion, color=Recur_status, fill=Recur_status)) +
  geom_boxplot(notch=FALSE, width=0.75, alpha=0.6, fill='white') +
  geom_quasirandom(dodge.width = 0.75, pch=21, col='black') +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
        legend.position = 'top') +
  labs(x='Cell Neighborhood', y='Proportion by Patient', 
       col='Recur Status', fill='Recur Status')
ggsave(file.path(OUTPUT_DIR, "CN_proportion_by_patient_boxplot_MSI-H_recur.pdf"),
       width = 7, height=6)

# MSS Recur vs No Recur ---------------------------------------------------

cn_props_long_msi_recur %>%
  filter(MSI == 'MSS') %>%
  ggplot(mapping = aes(x=CN, y=proportion, color=Recur_status, fill=Recur_status)) +
  geom_boxplot(notch=FALSE, width=0.75, alpha=0.6, fill='white') +
  geom_quasirandom(dodge.width = 0.75, pch=21, col='black') +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
        legend.position = 'top') +
  labs(x='Cell Neighborhood', y='Proportion by Patient', 
       col='Recur Status', fill='Recur Status')
ggsave(file.path(OUTPUT_DIR, "CN_proportion_by_patient_boxplot_MSS_recur.pdf"),
       width = 7, height=6)
