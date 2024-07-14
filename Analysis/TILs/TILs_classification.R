library(tidyverse)

set.seed(42)

sessionInfo()

# Load data ---------------------------------------------------------------

DATA_DIR <- './CRC_Spatial_Landscape/Data/'
OUTPUT_DIR <- DATA_DIR

cells <- read_csv(file.path(DATA_DIR, 'single_cell_metadata.csv'))
cn <- read_csv(file.path(DATA_DIR, 'cell_neighbourhoods.csv'))

# Preprocessing -----------------------------------------------------------

## Filter to cancer samples ---------------------------------------------

cells <- cells %>%
  filter(batch == 'Cancer')

## Merge single cell metadata and CNs --------------------------------------

roi_metadata <- cells %>% 
  group_by(ImageNumber) %>%
  filter(row_number()==1)

cn_metadata <- cn %>%
  inner_join(roi_metadata, by=c('File Name'='ImageNumber'))

# Classify lymphocytes and TILs -------------------------------------------

## Lymphocyte proportions --------------------------------------------------

# ROI level
lymphocytes_roi <- cn_metadata %>%
  select(alt_identifier, `File Name`, ClusterName) %>%
  group_by(alt_identifier, `File Name`, ClusterName) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = ClusterName, values_from = proportion) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  rowwise() %>%
  mutate(lymphocytes = sum(`CD4 T cell`, `CD8 T cell`, `B cell`)) %>%
  select(alt_identifier, `File Name`, lymphocytes)

# Patient level
lymphocytes_pt <- cn_metadata %>%
  select(alt_identifier, ClusterName) %>%
  group_by(alt_identifier, ClusterName) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = ClusterName, values_from = proportion) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  rowwise() %>%
  mutate(lymphocytes = sum(`CD4 T cell`, `CD8 T cell`, `B cell`)) %>%
  select(alt_identifier, lymphocytes)

## Intratumoural TILs proportion -------------------------------------------

# ROI level
itils_roi <- cn_metadata %>%
  filter(k8_cn %in% c('Bulk Tumor', 'P53+ Tumor', 'Proliferative Tumor')) %>%
  select(alt_identifier, `File Name`, ClusterName) %>%
  group_by(alt_identifier, `File Name`, ClusterName) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = ClusterName, values_from = proportion) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  rowwise() %>%
  mutate(itils = sum(`CD4 T cell`, `CD8 T cell`, `B cell`)) %>%
  select(alt_identifier, `File Name`, itils)
  
# Patient level
itils_pt <- cn_metadata %>%
  filter(k8_cn %in% c('Bulk Tumor', 'P53+ Tumor', 'Proliferative Tumor')) %>% 
  select(alt_identifier, ClusterName) %>%
  group_by(alt_identifier, ClusterName) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = ClusterName, values_from = proportion) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  rowwise() %>%
  mutate(itils = sum(`CD4 T cell`, `CD8 T cell`, `B cell`)) %>%
  select(alt_identifier, itils)

## Stromal TILs proportion -------------------------------------------------

# ROI level
stils_roi <- cn_metadata %>%
  filter(k8_cn %in% c('Immune-enriched stromal', 'Stromal')) %>%
  select(alt_identifier, `File Name`, ClusterName) %>%
  group_by(alt_identifier, `File Name`, ClusterName) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = ClusterName, values_from = proportion) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  rowwise() %>%
  mutate(stils = sum(`CD4 T cell`, `CD8 T cell`, `B cell`)) %>%
  select(alt_identifier, `File Name`, stils)

# Patient level
stils_pt <- cn_metadata %>%
  filter(k8_cn %in% c('Immune-enriched stromal', 'Stromal')) %>% 
  select(alt_identifier, ClusterName) %>%
  group_by(alt_identifier, ClusterName) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  select(-n) %>%
  pivot_wider(names_from = ClusterName, values_from = proportion) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  rowwise() %>%
  mutate(stils = sum(`CD4 T cell`, `CD8 T cell`, `B cell`)) %>%
  select(alt_identifier, stils)

# Merge tables ------------------------------------------------------------

# ROI level
tils_roi <- lymphocytes_roi %>%
  full_join(itils_roi, by=c('alt_identifier', 'File Name')) %>%
  full_join(stils_roi, by=c('alt_identifier', 'File Name')) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) 

# Patient level

tils_pt <- lymphocytes_pt %>%
  full_join(itils_pt, by='alt_identifier') %>%
  full_join(stils_pt, by='alt_identifier') %>%
  mutate_if(is.numeric, ~replace_na(., 0)) 

# Export tables -----------------------------------------------------------

tils_roi %>% 
  write_csv(file.path(OUTPUT_DIR, 'TILs_proportion_per_ROI.csv'))

tils_pt %>% 
  write_csv(file.path(OUTPUT_DIR, 'TILs_proportion_per_patient.csv'))