library(ggbeeswarm)
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

## Transform to long format ------------------------------------------------

cn_props_long <- cn_props %>%
  pivot_longer(`Bulk Tumor`:Stromal, names_to = 'CN', values_to = 'proportion')

## Order CNs by average abundance ------------------------------------------

cn_order_sorted <- cn_props_long %>% 
  group_by(CN) %>%
  summarise(median = median(proportion)) %>%
  arrange(median) %>%
  pull(CN)

cn_props_long <- cn_props_long %>%
  mutate(CN = factor(CN, levels = cn_order_sorted))

# Boxplot -----------------------------------------------------------------

cn_props_long %>%
  ggplot(mapping = aes(x=proportion, y=CN, col=CN, fill=CN)) + 
  geom_quasirandom(pch=21, col='gray10') +
  geom_boxplot(alpha=0) + 
  labs(x = 'Proportion of CN in each patient', y='') +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values = c('#A6D854', '#E5C494', '#66C2A5', '#FC8D62', '#D62728', '#FFD92F', '#E78AC3', '#8DA0CB'),
                     breaks = c('P53+ Tumor', 'Stromal', 'Bulk Tumor', 
                                'Immune-enriched stromal', 'CD8+ T cell enriched', 
                                'Proliferative Tumor', 'NK-cell enriched', 'Mixed Immune')) +
  scale_fill_manual(values = c('#A6D854', '#E5C494', '#66C2A5', '#FC8D62', '#D62728', '#FFD92F', '#E78AC3', '#8DA0CB'),
                     breaks = c('P53+ Tumor', 'Stromal', 'Bulk Tumor', 
                                'Immune-enriched stromal', 'CD8+ T cell enriched', 
                                'Proliferative Tumor', 'NK-cell enriched', 'Mixed Immune')) +
  theme_bw() +
  theme(legend.position = 'none')

ggsave(file.path(OUTPUT_DIR, 'CN_proportion_by_patient_summary_boxplot_beeswarm.pdf'), width=7, height=4)
