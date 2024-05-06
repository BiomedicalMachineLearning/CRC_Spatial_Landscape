# library(boot)
# library(circlize)
# library(corrplot)
# library(data.table)
# library(Hmisc)
library(ggbeeswarm)
# library(limma)
# library(RColorBrewer)
# library(rstatix)
# library(Seurat)
# library(speckle)
library(readxl)
library(tidyverse)

set.seed(42)

sessionInfo()


# Load data ---------------------------------------------------------------

TCGA_DIR <- './CRC_Spatial_Landscape/Analysis/TCGA_validation'

CIBERSORT_PATH <- file.path(TCGA_DIR, 'CIBERSORTx_Results.txt')
METADATA_PATH <- file.path(TCGA_DIR, 'Clinical_metadata.xlsx')
MSI_PATH <- file.path(TCGA_DIR, 'MSI_cBioPortal.tsv')

OUTPUT_DIR <- TCGA_DIR

cibersort <- read_tsv(CIBERSORT_PATH)
msi <- read_tsv(MSI_PATH) 
meta <- read_excel(METADATA_PATH, 
                   sheet = 'TCGA-CDR',
                   na = c('', '[Not Available]', '[Not Applicable]'))

# Filter to stage III colon adenocarcinoma --------------------------------

## Filter metadata ---------------------------------------------------------

# Filter to only colon adenocarcinomas
meta_coad <- meta %>%
  filter(type == 'COAD')

n_distinct(meta_coad$bcr_patient_barcode)

unique(meta_coad$ajcc_pathologic_tumor_stage)

# Filter to only stage III cancers
meta_coad_III <- meta_coad %>%
  filter(grepl('III', ajcc_pathologic_tumor_stage)) %>%
  separate(bcr_patient_barcode, into = c('TCGA', 'TSS', 'patient'), remove = FALSE) %>%
  unite(col = 'TSS_patient', TSS:patient, sep='-')

unique(meta_coad_III$ajcc_pathologic_tumor_stage)

meta_coad_III

nrow(meta_coad_III)

n_distinct(meta_coad_III$bcr_patient_barcode)

## Filter CIBERSORTx -------------------------------------------------------

### Filter to cancer samples -----------------------------------------------

# View(cibersort)

cibersort <- cibersort %>% 
  separate(Mixture, 
           into = c('TCGA', 'TSS', 'patient', 'sample', 'portion', 'plate', 'center'),
           remove = FALSE) %>%
  unite(col = 'TSS_patient', TSS:patient, sep='-')

table(cibersort$sample, useNA = 'ifany')

# Tumor sample codes begin with a 0 whereas normal begins with 1
cibersort_ca <- cibersort %>%
  filter(grepl("0\\d\\w", sample))

table(cibersort_ca$sample, useNA = 'ifany')


### Merge duplicate patients ------------------------------------------------

dup_pts <- (table(cibersort_ca$TSS_patient) > 1)%>%
  as_tibble(rownames = 'TSS_patient') %>%
  filter(value == TRUE) %>%
  pull(TSS_patient)

n_distinct(dup_pts)

cibersort_ca %>% 
  filter(TSS_patient %in% dup_pts) %>%
  arrange(TSS_patient)

table(cibersort_ca$P.value, useNA = 'ifany')
table(cibersort_ca$Correlation, useNA = 'ifany')
table(cibersort_ca$RMSE, useNA = 'ifany')

cibersort_ca_long <- cibersort_ca %>% 
  pivot_longer(cols = Fibroblasts:PMNs,
               names_to = 'cell_type', 
               values_to = 'proportion') %>%
  select(Mixture, TSS_patient, cell_type, proportion)

# Check that all proportions add up to 1
cibersort_ca_long %>%
  group_by(Mixture) %>%
  summarise(total_prop = sum(proportion)) %>%
  filter(total_prop < 0.99 | total_prop > 1.01)

# Where more than one sample exists for a patient, replace with mean of the samples' cell type proportions
cibersort_ca_long_mean <- cibersort_ca_long %>%
  select(-Mixture) %>%
  group_by(TSS_patient, cell_type) %>%
  summarise(proportion = mean(proportion))

# Check that all proportions add up to 1
cibersort_ca_long_mean %>% 
  group_by(TSS_patient) %>%
  summarise(total_prop = sum(proportion)) %>%
  filter(total_prop < 0.99 | total_prop > 1.01)

cibersort_ca_wide <- cibersort_ca_long_mean %>%
  pivot_wider(names_from = cell_type,
              values_from = proportion)

cibersort_ca_wide

### Filter to stage III COAD ------------------------------------------------

cibersort_coad_III <- cibersort_ca_wide %>%
  inner_join(meta_coad_III %>% select(-TCGA),
             by = 'TSS_patient')

cibersort_coad_III

nrow(cibersort_coad_III)

n_distinct(cibersort_coad_III$TSS_patient)


## MSI determination -------------------------------------------------------

msi_wide <- msi %>% 
  pivot_longer(cols = contains('TCGA-'),
               names_to = 'TCGA_barcode') %>%
  pivot_wider(names_from = track_name, 
              values_from = value) %>%
  drop_na(`MSIsensor Score`) %>%
  mutate(MSIsensor = as.double(`MSIsensor Score`)) %>%
  arrange(MSIsensor) %>%
  mutate(TCGA_barcode = factor(TCGA_barcode, levels = TCGA_barcode)) %>%
  mutate(MSI_status = ifelse(MSIsensor > 10, 'MSI-H', 'MSS'),
         MSI_status = factor(MSI_status, levels = c('MSI-H', 'MSS'))) %>%
  separate(TCGA_barcode, 
           into=c('TCGA', 'TSS', 'patient'),
           remove = FALSE) %>%
  unite(col = 'TSS_patient', TSS:patient, sep='-')

nrow(msi_wide) 

any(table(msi_wide$TCGA_barcode) > 1)


unique(msi_wide$Subtype)

msi_wide %>%
  ggplot(aes(y=MSIsensor, x=TCGA_barcode, col=MSI_status)) +
  geom_point(size=0.5) +
  geom_hline(yintercept = 10, linetype=2) +
  theme(axis.text.x = element_blank()) + 
  labs(x='TCGA Patients', y='MSIsensor score',
       subtitle = 'TCGA-COAD cohort')
ggsave(file.path(OUTPUT_DIR, 'MSIsensor_scatter_plot.pdf'), height = 3, width = 6)

msi_wide

## Merge CIBERSORT, metadata and MSI -------------------------------------

coad_III <- cibersort_coad_III %>%
  left_join(msi_wide, by='TSS_patient') %>%
  mutate(Recur_status = case_match(DFI,
                                   0 ~ 'no recur',
                                   1 ~ 'recur', 
                                   NA ~ NA),
         Recur_status = factor(Recur_status, levels = unique(Recur_status))) %>%
  ungroup()

nrow(coad_III)

table(coad_III$MSI_status, useNA = 'ifany')

unique(coad_III$Recur_status)


# Cell type proportions ---------------------------------------------------

coad_III_long <- coad_III %>%
  pivot_longer(cols = B.cells:PMNs,
               names_to = 'cell_type',
               values_to = 'proportion')

## Plot overall cell type proportion ---------------------------------------

ct_order <- coad_III_long %>%
  group_by(cell_type) %>%
  summarise(median_props = median(proportion)) %>%
  arrange(median_props) %>%
  pull(cell_type)

coad_III_long$cell_type <- factor(coad_III_long$cell_type,
                                  levels = ct_order)

ct_cols <- brewer.pal(12, 'Set3')
names(ct_cols) <- ct_order
ct_cols

# All cell types
coad_III_long %>%
  ggplot(aes(x=proportion, y=cell_type, col=cell_type, fill=cell_type)) +
  geom_boxplot(alpha = 0.6, fill='white') +
  geom_quasirandom(pch=21, col='black', size=0.8, stroke=0.2) +
  labs(x='Proportion of total cells',
       y='Cell types') +
  theme_bw() +
  theme(axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black'),
        legend.position = 'none') +
  scale_color_manual(values = ct_cols) +
  scale_fill_manual(values = ct_cols)
ggsave(file.path(OUTPUT_DIR, 'TCGA_cell_type_proportions_boxplot.pdf'),
       width=5, height=4)

# Immune cell types
coad_III_long %>%
  filter(!cell_type %in% c('Epithelial.cells', 'Fibroblasts', 'Endothelial.cells')) %>%
  ggplot(aes(x=proportion, y=cell_type, col=cell_type, fill=cell_type)) +
  geom_boxplot(alpha = 0.6, fill='white') +
  geom_quasirandom(pch=21, col='black', size=0.8, stroke=0.2) +
  labs(x='Proportion of total cells',
       y='Immune cell types') +
  theme_bw() +
  theme(axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black'),
        legend.position = 'none') +
  scale_color_manual(values = ct_cols) +
  scale_fill_manual(values = ct_cols)
ggsave(file.path(OUTPUT_DIR, 'TCGA_immune_cell_type_proportions_boxplot.pdf'),
       width=5, height=4)




# Export ------------------------------------------------------------------

coad_III %>%
  write_csv(file.path(OUTPUT_DIR, 'TCGA_COAD_Stage_III_CIBERSORTx.csv'))