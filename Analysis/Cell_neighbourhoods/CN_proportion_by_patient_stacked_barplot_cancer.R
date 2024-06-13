library(circlize)
library(ComplexHeatmap)
library(tidyverse)

sessionInfo()

set.seed(42)

# Load data ---------------------------------------------------------------

DATA_DIR <- './CRC_Spatial_Landscape/Data/'
OUTPUT_DIR <- './CRC_Spatial_Landscape/Analysis/Cell_neighbourhoods'

samples <- read_csv(file.path(DATA_DIR, 'samples.csv'))

cn_props <- read_csv(file.path(DATA_DIR, "cell_neighbourhoods_k8_proportions_by_patient.csv"))


# Plotting parameters -----------------------------------------------------

cn_order = c('CD8+ T cell enriched', 
             'Mixed Immune',
             'NK-cell enriched',
             'Immune-enriched stromal',
             'Stromal',
             'P53+ Tumor',
             'Proliferative Tumor',
             'Bulk Tumor')

cn_colours <-  c('Bulk Tumor' = 'yellow',
                 'Proliferative Tumor' ='red',
                 'P53+ Tumor' = 'orange',
                 'Immune-enriched stromal' = 'blue',
                 'CD4 T cell' = 'light green',
                 'CD8+ T cell enriched' = 'dark green',
                 'Macrophages' = 'light blue',
                 'NK-cell enriched' = 'purple',
                 'Mixed Immune' = 'pink',
                 'Stromal' = 'grey')

# Preprocess cell neighbourhood proportions -------------------------------

## Filter for cancer samples -----------------------------------------------

cn_props <- cn_props %>%
  filter(batch=='Cancer')

## Order patients by abundance of CD8+ T cell enriched CNs -----------------

cn_props <- cn_props %>%
  arrange(desc(`CD8+ T cell enriched`))


## Create patient x CN matrix ----------------------------------------------

cn_props_mat <- cn_props %>%
  select(all_of(cn_order)) %>%
  as.matrix()

rownames(cn_props_mat) <- cn_props %>% pull(`alt_identifier`)

# Plot stacked barplot ----------------------------------------------------

cols <- cn_colours[cn_order]

# Plot barplot as a heatmap annotation
cn_bar <- HeatmapAnnotation(`CN Proportions` = anno_barplot(cn_props_mat,
                                               gp = gpar(fill = cols),
                                               bar_width = 1,
                                               height = unit(6, 'cm')))

cn_legend <-Legend(labels = colnames(cn_props_mat),
                   legend_gp = gpar(fill = cols),
                   title = 'Cell Neighbourhoods',
                   nrow = 2)

# Create a placeholder heatmap of zero height
cn_ht <-  Heatmap(t(cn_props_mat), 
                  name ='mat', 
                  rect_gp = gpar(type = 'none'),
                  cluster_columns = FALSE,
                  show_column_dend=FALSE,
                  show_row_names = FALSE,
                  show_row_dend = FALSE,
                  column_title = 'Patients',
                  column_title_side = 'bottom',
                  height = unit(0, 'cm'),
                  show_heatmap_legend = FALSE,
                  top_annotation= cn_bar,
                  column_names_gp = gpar(fontsize=10))

## Export stacked barplot --------------------------------------------------

pdf(file.path(OUTPUT_DIR, "CN_proportion_stacked_barplot_cancer.pdf"), width = 10, height = 5)
draw(cn_ht, annotation_legend_list = list(cn_legend), annotation_legend_side = 'bottom')
dev.off()

