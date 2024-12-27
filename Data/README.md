[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13901180.svg)](https://doi.org/10.5281/zenodo.13901180)

# Data Overview

We acquired imaging mass cytometry (IMC) data for 52 patients with stage III colon cancer. The raw image files are availabe on our [Zenodo](https://doi.org/10.5281/zenodo.13901180). We transformed these multiplexed images into single cell protein expression data whilst maintaining spatial coordinates. On our [Zenodo](https://doi.org/10.5281/zenodo.13901180) these are provided are Anndata (`.h5ad`) and Seurat (`.rds`) objects for analysis in Python and R respectively. 

Integrated with the IMC data are matched Whole-Exome-Sequencing (WES) data for microsatellite instability (MSI) determination and H&E stained histopathological whole-slide-images (WSIs) for pathologist validation of cell types. We provide additional metadata such as patient demographics and clinical outcomes (i.e. overall survival, recurrence) in the `samples.csv` file. 

We validated our findings against [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/ccg/research/genome-sequencing/tcga) and from CODEX spatial proteomics data of colorectal cancer by [Schurch et. al.](https://doi.org/10.1016/j.cell.2020.07.005)

# Raw data
## Sample metadata

`samples.csv` contains the metadata for each sample including:
- MSI-status
- Clinical outcomes (e.g. survival, recurrence)

It also contains the cell type proportion for each of the 10 cell types.

# Processed data

## Single cell objects

Seurat and Anndata objects are released on our [Zenodo](https://doi.org/10.5281/zenodo.13901180). In this repo we provide views of this data so that our scripts can be more easily run. This includes `single_cell_phenotypes_and_coordinates.csv`. 

# Analysis outputs

Although analysis outputs can be produced by running the code, we provide relevant outputs for ease of downstream analysis. 

## Cell type proportions

We summarise the proportional abundance of each cell type relative to the total no. of cells per region-of-interest (ROI) or per patient in `CRC_Spatial_Landscape/Data/cell_type_proportion_per_CN_by_ROI.csv` and `CRC_Spatial_Landscape/Data/cell_type_proportion_per_CN_by_patient.csv` respectively.

## Cell neighborhoods

`cell_neighbourhoods_k8_proportions_by_patient.csv` contains the proportion of cell neighbourhoods for each patient.

## Tumor-infiltrating lymphocytes (TILs)

ROI level: `CRC_Spatial_Landscape/Data/TILs_proportion_per_ROI.csv`
Patient level: `CRC_Spatial_Landscape/Data/TILs_proportion_per_patient.csv`
