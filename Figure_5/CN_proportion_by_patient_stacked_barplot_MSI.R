library(circlize)
library(ComplexHeatmap)
library(tidyverse)

sessionInfo()

set.seed(42)

# Load data ---------------------------------------------------------------

DATA_DIR <- './CRC_Spatial_Landscape/Data/'
OUTPUT_DIR <- './CRC_Spatial_Landscape/Figure_5'

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
  filter(batch=='Cancer') %>%
  left_join(samples %>% select(alt_identifier,
                               MSI_any,
                               Reccurrance),
            by='alt_identifier')

## Order patients by abundance of CD8+ T cell enriched CNs -----------------

cn_props <- cn_props %>%
  arrange(desc(`CD8+ T cell enriched`))


## Create patient x CN matrix for MSI and MSS samples ----------------------

### MSI-H ----------------------------------------------------------------

cn_props_msi_pos <- cn_props %>% filter(MSI_any == TRUE)

cn_props_msi_pos_mat <- cn_props_msi_pos %>%
  select(all_of(cn_order)) %>%
  as.matrix()

rownames(cn_props_msi_pos_mat) <- cn_props_msi_pos %>% pull(alt_identifier)

### MSS ------------------------------------------------------------------

cn_props_msi_neg <- cn_props %>% filter(MSI_any == FALSE)

cn_props_msi_neg_mat <- cn_props_msi_neg %>%
  select(all_of(cn_order)) %>%
  as.matrix()

rownames(cn_props_msi_neg_mat) <- cn_props_msi_neg %>% pull(alt_identifier)


# Plot stacked barplots ---------------------------------------------------

cols <- cn_colours[cn_order]


## MSI-H ------------------------------------------------------------------

# Annotate which patients recurred or not
col_ha_msi_pos <- HeatmapAnnotation(df = as.data.frame(cn_props_msi_pos %>% select(Reccurrance)),
                                col = list(Reccurrance = c('no recur'='green', 'recur'='black',
                                                           'unk'='gray', 'unk - no FU'='gray',
                                                           'unk, moved'='gray')))

# Plot barplot as a heatmap annotation
stacked_bar_msi_pos <-  HeatmapAnnotation(`MSI-H` = anno_barplot(cn_props_msi_pos_mat, 
                                                                 gp = gpar(fill = cols),
                                                                 bar_width = 1, 
                                                                 height = unit(10, "cm")),
                                          show_legend=c(TRUE, TRUE))

lgd_msi_pos <- Legend(labels = colnames(cn_props_msi_pos_mat), 
              legend_gp = gpar(fill = cols), 
              title = 'Cell Neighborhoods',
              nrow = 2)

# Create a placeholder heatmap of zero height
ht_msi_pos <-  Heatmap(t(cn_props_msi_pos_mat), 
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
               top_annotation = stacked_bar_msi_pos,
               column_names_gp = gpar(fontsize=10))

### Export stacked barplot -------------------------------------------------
pdf(file.path(OUTPUT_DIR, "CN_proportion_stacked_barplot_no_legend_MSI.pdf"), width = 4, height = 8)
draw(col_ha_msi_pos %v% ht_msi_pos)
dev.off()

## MSS --------------------------------------------------------------------

# Annotate which patients recurred or not
col_ha_msi_neg <- HeatmapAnnotation(df = as.data.frame(cn_props_msi_neg %>% select(Reccurrance)),
                                    col = list(Reccurrance = c('no recur'='green', 'recur'='black',
                                                               'unk'='gray', 'unk - no FU'='gray',
                                                               'unk, moved'='gray')))

# Plot barplot as a heatmap annotation
stacked_bar_msi_neg <-  HeatmapAnnotation(MSS = anno_barplot(cn_props_msi_neg_mat, 
                                                                 gp = gpar(fill = cols),
                                                                 bar_width = 1, 
                                                                 height = unit(10, "cm")),
                                          show_legend=c(TRUE, TRUE))

lgd_msi_neg <- Legend(labels = colnames(cn_props_msi_neg_mat), 
                      legend_gp = gpar(fill = cols), 
                      title = 'Cell Neighborhoods',
                      nrow = 2)

# Create a placeholder heatmap of zero height
ht_msi_neg <-  Heatmap(t(cn_props_msi_neg_mat), 
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
                       top_annotation = stacked_bar_msi_neg,
                       column_names_gp = gpar(fontsize=10))

### Export stacked barplot -------------------------------------------------
pdf(file.path(OUTPUT_DIR, "CN_proportion_stacked_barplot_no_legend_MSS.pdf"), width = 8, height = 8)
draw(col_ha_msi_neg %v% ht_msi_neg)
dev.off()
