library(ggbeeswarm)
library(ggpubr)
library(tidyverse)

set.seed(42)

sessionInfo()


# Load data ---------------------------------------------------------------

DATA_DIR <- './CRC_Spatial_Landscape/Data/'
OUTPUT_DIR <- './CRC_Spatial_Landscape/Analysis/Dataset_overview'

cells <- read_csv(file.path(DATA_DIR, 'single_cell_metadata.csv'))


# Filter to cancer samples ------------------------------------------------

cells <- cells %>%
  filter(batch == 'Cancer')

# Compute tumour purity ---------------------------------------------------

## Tumour purity per ROI ---------------------------------------------------

# Calculate proportion of each cell type per ROI
props_roi <- cells %>%
  select(alt_identifier, ImageNumber, cell_type) %>%
  group_by(alt_identifier, ImageNumber, cell_type) %>%
  summarise(n= n()) %>%
  mutate(proportion = n / sum(n))

# Calculate tumour purity by summing proportion of tumour cell types
props_roi_wide <- props_roi %>%
  select(-n) %>%
  pivot_wider(names_from = cell_type, values_from = proportion) %>%
  rowwise() %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  mutate(tumour_purity = sum(`Epithelial/Tumor`,
                             `P53+ Tumor`, 
                             `Proliferative Tumor`))

## Tumour purity per patient -----------------------------------------------

# Calculate proportion of each cell type per patient
props_pt <- cells %>%
  select(alt_identifier, cell_type) %>%
  group_by(alt_identifier, cell_type) %>%
  summarise(n= n()) %>%
  mutate(proportion = n / sum(n))

# Calculate tumour purity by summing proportion of tumour cell types

props_pt_wide <- props_pt %>%
  select(-n) %>%
  pivot_wider(names_from = cell_type, values_from = proportion) %>%
  rowwise() %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  mutate(tumour_purity = sum(`Epithelial/Tumor`,
                             `P53+ Tumor`, 
                             `Proliferative Tumor`))

# Plot tumour purity distributions ----------------------------------------

purity_boxplot_roi <- props_roi_wide %>%
  mutate(x='') %>%
  ggplot(aes(x=x, y=tumour_purity)) +
  geom_boxplot(outlier.shape=NA) +
  geom_quasirandom(size=0.8) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x='ROI level',
       y='Tumor Purity') +
  ylim(0,1)

purity_boxplot_pt <- props_pt_wide %>%
  mutate(x='') %>%
  ggplot(aes(x=x, y=tumour_purity)) +
  geom_boxplot(outlier.shape=NA) +
  geom_quasirandom(size=0.8) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x='Patient level',
       y='Tumor Purity') +
  ylim(0,1)

ggarrange(purity_boxplot_roi, purity_boxplot_pt + rremove('ylab'))
ggsave(file.path(OUTPUT_DIR, 'tumor_purity_boxplots_roi_vs_patient.pdf'),
       width = 4, height=4)


# Plot tumour purity distributions within each patient --------------------

# Order patients by median tumour purity
alt_identifier_order <- props_roi_wide %>%
  group_by(alt_identifier) %>%
  summarise(median_tumour_purity = median(tumour_purity)) %>%
  arrange(median_tumour_purity) %>%
  pull(alt_identifier)

props_roi_wide %>%
  mutate(alt_identifier = factor(alt_identifier, levels=alt_identifier_order)) %>%
  ggplot(aes(x=alt_identifier, y=tumour_purity)) +
  geom_boxplot(width=0.75, alpha=0.6, fill='white') +
  geom_quasirandom(dodge.width=0.75, pch=21, col='black') +
  scale_x_discrete(guide = guide_axis(angle=45)) +
  theme_bw() +
  labs(x = 'Patient ID',
       y = 'Tumor purity per ROI')

ggsave(file.path(OUTPUT_DIR, 'Tumor_purity_per_ROI_by_patient_boxplot.pdf'), width=9, height=4.5)

# Export tumour purity tables ---------------------------------------------

props_roi_wide %>%
  write_csv(file.path(OUTPUT_DIR, 'ROI_tumor_purity.csv'))

props_pt_wide %>%
  write_csv(file.path(OUTPUT_DIR, 'patient_tumor_purity.csv'))
