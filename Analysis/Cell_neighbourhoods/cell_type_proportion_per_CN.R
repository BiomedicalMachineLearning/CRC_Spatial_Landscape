library(tidyverse)

set.seed(42)

sessionInfo()

# Load data ---------------------------------------------------------------

DATA_DIR <- './CRC_Spatial_Landscape/Data/'
OUTPUT_DIR <- DATA_DIR

cells <- read_csv(file.path(DATA_DIR, 'single_cell_metadata.csv'))
cn <- read_csv(file.path(DATA_DIR, 'cell_neighbourhoods.csv'))

# Merge single cell metadata and CNs --------------------------------------

roi_metadata <- cells %>% 
  group_by(ImageNumber) %>%
  filter(row_number()==1)

cn_metadata <- cn %>%
  inner_join(roi_metadata, by=c('File Name'='ImageNumber'))

# Calculate cell type proportions per CN ----------------------------------

## ROI level ---------------------------------------------------------------

props_roi <- cn_metadata %>%
  select(alt_identifier, `File Name`, k8_cn, ClusterName) %>%
  group_by(alt_identifier, `File Name`, k8_cn, ClusterName) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n))

# Transform into wide format
props_roi_wide <- props_roi %>%
  select(-n) %>%
  pivot_wider(names_from = ClusterName, values_from = proportion) %>%
  mutate_if(is.numeric, ~replace_na(., 0))

## Patient level --------------------------------------------------------

props_pt <- cn_metadata %>%
  select(alt_identifier, k8_cn, ClusterName) %>%
  group_by(alt_identifier, k8_cn, ClusterName) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n))

# Transform into wide format
props_pt_wide <- props_pt %>%
  select(-n) %>%
  pivot_wider(names_from = ClusterName, values_from = proportion) %>%
  mutate_if(is.numeric, ~replace_na(., 0))

# Export tables -----------------------------------------------------------

props_roi_wide %>% 
  write_csv(file.path(OUTPUT_DIR, 'cell_type_proportion_per_CN_by_ROI.csv'))

props_pt_wide %>% 
  write_csv(file.path(OUTPUT_DIR, 'cell_type_proportion_per_CN_by_patient.csv'))