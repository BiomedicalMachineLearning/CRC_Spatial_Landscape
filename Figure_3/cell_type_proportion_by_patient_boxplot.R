library(circlize)
library(ggbeeswarm)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)

sessionInfo()

set.seed(42)

# Load cell type proportions ----------------------------------------------

SAMPLES_PATH <- './CRC_Spatial_Landscape/Data/samples.csv'
OUTPUT_DIR <- './CRC_Spatial_Landscape/Figure_3'

samples <- read_csv(SAMPLES_PATH)

# Plotting parameters -----------------------------------------------------

ct_colours <-  c('Epithelial/Tumor' = 'yellow',
                 'Proliferative Tumor' ='red',
                 'P53+ Tumor' = 'orange',
                 'B cell' = 'blue',
                 'CD4 T cell' = 'light green',
                 'CD8 T cell' = 'dark green',
                 'Macrophages' = 'light blue',
                 'NK-cells' = 'purple',
                 'Proliferative Immune' = 'pink',
                 'Stromal' = 'grey')

# Boxplots ----------------------------------------------------------------

## All cell types --------------------------------------------------------

all_cell_types <- samples %>% 
  pivot_longer(cols = `B cell`:Stromal, names_to = 'cell_type', values_to = 'proportion') %>%
  mutate(cell_type = as.factor(cell_type))

# Order cell types by abundance
ct_order_sorted <- all_cell_types %>% 
  group_by(cell_type) %>%
  summarise(median = median(proportion)) %>%
  arrange(median) %>%
  pull(cell_type)

all_cell_types <- all_cell_types %>% 
  mutate(cell_type = factor(cell_type, levels = ct_order_sorted))

boxplot_all <- all_cell_types %>%
  ggplot(aes(x=proportion, y=cell_type, col=cell_type, fill=cell_type)) +
  geom_boxplot(alpha=0.6, fill='white') +
  geom_quasirandom(pch=21, col='black') +
  labs(x = 'Proportion of total cells') +
  theme(axis.title.y = element_blank(),
        legend.position = 'none') +
  scale_color_manual(values = ct_colours) +
  scale_fill_manual(values = ct_colours)

## Immune cell types -----------------------------------------------------

immune_subtypes <- c('CD8 T cell', 'CD4 T cell', 'B cell', 'Macrophages', 
                     'NK-cells', 'Proliferative Immune')

# Scale immune cell type proportions relative to all immune cells in a patient sample
immune_mat <- samples %>%
  select(all_of(immune_subtypes)) %>%
  as.matrix() 

immune_mat_scaled <- immune_mat * (1 / rowSums(immune_mat))

immune_cell_types <- as_tibble(immune_mat_scaled, rownames='alt_identifier') %>%
  pivot_longer(cols=`CD8 T cell`:`Proliferative Immune`, names_to='cell_type', values_to = 'proportion')

# Order cell types by adundance
immune_ct_order_sorted <- immune_cell_types %>% 
  group_by(cell_type) %>%
  summarise(median = median(proportion)) %>%
  arrange(median) %>%
  pull(cell_type)

immune_cell_types <- immune_cell_types %>% 
  mutate(cell_type = factor(cell_type, levels = immune_ct_order_sorted))

boxplot_immune <- immune_cell_types %>%
  ggplot(aes(x=proportion, y=cell_type, col=cell_type, fill=cell_type)) +
  geom_boxplot(alpha=0.6, fill='white') +
  geom_quasirandom(pch=21, col='black') +
  labs(x = 'Proportion of immune cells') +
  theme(axis.title.y = element_blank(),
        legend.position = 'none') +
  scale_color_manual(values = ct_colours) +
  scale_fill_manual(values = ct_colours)


## Export combined boxplot -----------------------------------------------

ggarrange(boxplot_all, boxplot_immune, ncol = 2, nrow = 1)
ggsave(file.path(OUTPUT_DIR, 'cell_type_proportion_by_patient_boxplot_total_and_immune.pdf'), width = 7, height = 5)


