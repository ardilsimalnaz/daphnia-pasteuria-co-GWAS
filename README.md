# Daphnia–Pasteuria co-GWAS

Scripts used for population structure and co-genome-wide association analyses of
*Daphnia magna* and *Pasteuria ramosa* as part of my bachelor thesis at the
Technical University of Munich (TUM).

## Repository structure

### 01_population_structure

Scripts used to investigate population structure in *D. magna* and *P. ramosa*.

- `PCA/`
  - `plot_daphnia_pca.R`
  - `plot_pasteuria_pca.R`

- `ADMIXTURE/`
  - `plot_daphnia_admixture_K6.R`
  - `plot_pasteuria_admixture_K3.R`
  - `plot_pasteuria_admixture_K9_final.R`
  - `plot_admixture_deltaCV.R`

- `fastSTRUCTURE/`
  - `plot_daphnia_faststructure_K6.R`
  - `plot_pasteuria_faststructure_K3.R`
  - `structure.py`
  - `chooseK.py`

### 02_model1_cogwas

Scripts for the Model 1 co-GWAS workflow.

- `model1 full reproducible.sh`  
  Main reproducible Model 1 pipeline including genotype filtering, PCA,
  construction of binary *P. ramosa* phenotypes, and co-GWAS analysis.

- `cogwas_onevariant_pipeline.sh`  
  One-variant workflow used to validate the Model 1 analysis procedure.

### 03_figures

Scripts used to generate the main Model 1 co-GWAS visualizations.

- `plot_model1_manhattan_final.R`
- `plot_model1_circle_paperstyle.R`

### 04_resistance_loci_PCL

Post-hoc analyses of Model 1 associations involving known *D. magna* resistance
loci and *P. ramosa* candidate loci (PCLs).

- `model1_resistance_analysis.sh`  
  Identifies Model 1 associations overlapping the ABC/C, D and E resistance
  loci and tests the corresponding *P. ramosa* variants for overlap with PCL
  triplets.

- `PCL_gene_LD_check.sh`  
  Assigns PCL-overlapping variants to individual PCL genes and calculates
  pairwise unphased genotype correlations (r²).

## Software

The analyses were performed using:

- PLINK 2.0
- R
- ADMIXTURE
- fastSTRUCTURE
- Bash / command-line utilities

## Data

Raw genotype files and large intermediate analysis files are not included in
this repository because of their size. The scripts therefore require the
corresponding filtered genotype and supporting input files described in the
bachelor thesis.

## Thesis

These scripts accompany the bachelor thesis on joint genomic analysis of the
host–parasite system *Daphnia magna*–*Pasteuria ramosa*.

Author: Simal Naz Ardil  
Technical University of Munich (TUM), 2026
