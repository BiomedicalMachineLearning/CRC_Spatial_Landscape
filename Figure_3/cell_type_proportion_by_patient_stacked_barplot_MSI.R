library(circlize)
library(ComplexHeatmap)
library(tidyverse)

sessionInfo()

set.seed(42)

# Load data ---------------------------------------------------------------

DATA_DIR <- './CRC_Spatial_Landscape/Data/'
OUTPUT_DIR <- './CRC_Spatial_Landscape/Figure_3'

samples <- read_csv(file.path(DATA_DIR, 'samples.csv'))

# Plotting parameters -----------------------------------------------------

cell_type_order <- c('CD8 T cell',
                     'CD4 T cell',
                     'B cell', 
                     'Macrophages',
                     'NK-cells', 
                     'Proliferative Immune', 
                     'Stromal', 
                     'P53+ Tumor', 
                     'Proliferative Tumor', 
                     'Epithelial/Tumor')

cell_type_colours <- c('Epithelial/Tumor' = 'yellow',
                       'Proliferative Tumor' ='red',
                       'P53+ Tumor' = 'orange',
                       'B cell' = 'blue',
                       'CD4 T cell' = 'light green',
                       'CD8 T cell' = 'dark green',
                       'Macrophages' = 'light blue',
                       'NK-cells' = 'purple',
                       'Proliferative Immune' = 'pink',
                       'Stromal' = 'grey')

cols <- cell_type_colours[cell_type_order]

# Preprocess cell type proportions and metadata ---------------------------

## Order patients by MSI status and CD8+ T cell abundance ----
samples <- samples %>%
  drop_na(MSI_any) %>%
  arrange(desc(MSI_any), desc(`CD8 T cell`))

## Create patient x cell types matrix ----

cell_type_mat <- samples %>% 
  select(all_of(cell_type_order)) %>%
  as.matrix()

rownames(cell_type_mat) <- samples %>% pull(alt_identifier)

## Consolidate recur status labels

samples <- samples %>%
  mutate(Recur_status = case_when(Reccurrance == 'recur' ~ 'Recur',
                                  Reccurrance == 'no recur' ~ 'No recur',
                                  str_detect(Reccurrance, 'unk') ~ 'Unknown'))

## Select metadata for plotting ----

metadata_df <- samples %>%
  select(CHR_17P_DEL, Recur_status, `5YearSurvival`) %>%
  as.data.frame()

# Create metadata annotations ---------------------------------------------

metadata_legend_params <- list(
  CHR_17P_DEL = list(
    at = c('YES', 'NO'),
    labels = c('Yes', 'No')
  ),
  `5YearSurvival` = list(
    at = c('<5', '>5'),
    labels = c('< 5 years', '> 5 years')
  )
)

metadata_ha <- HeatmapAnnotation(df = metadata_df,
                                 col = list(CHR_17P_DEL = c('YES' = 'red', 
                                                            'NO' = 'blue'),
                                            Recur_status = c('No recur' = 'green', 
                                                            'Recur' = 'black',
                                                            'Unknown' = 'gray'),
                                            `5YearSurvival` = c('<5' = 'black', 
                                                                '>5' = 'green')),
                                 annotation_label=c('Chr 17p del', 
                                                    'Recur status', 
                                                    '5 year survival'),
                                 annotation_legend_param = metadata_legend_params)

# Plot stacked barplot ----------------------------------------------------

# Plot barplot as a heatmap annotation
proportions_bar <- HeatmapAnnotation(Proportion = anno_barplot(cell_type_mat,
                                                               gp = gpar(fill = cols),
                                                               bar_width = 1, 
                                                               height = unit(10, 'cm')))

lgd <- Legend(labels = colnames(cell_type_mat),
              title = 'Cell Types',
              legend_gp = gpar(fill = cols),
              nrow = 2)

# Create a placeholder heatmap of zero height
ht <-  Heatmap(t(cell_type_mat), 
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
               top_annotation= proportions_bar,
               column_names_gp = gpar(fontsize=10),
               column_split = factor(samples$MSI_any, levels=c('TRUE', 'FALSE')),
               column_gap = unit(5, 'mm')
               )

## Export stacked barplot --------------------------------------------------

pdf(file.path(OUTPUT_DIR, "cell_type_proportion_stacked_barplot_MSI.pdf"), width = 10, height = 7)
draw(metadata_ha %v% ht, 
     annotation_legend_list = list(lgd), 
     annotation_legend_side = 'bottom')
dev.off()
