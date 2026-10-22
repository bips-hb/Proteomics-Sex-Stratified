# Proteomics Sex Stratified

R code for the sex-stratified proteomic analyses conducted in a subsample of the German National Cohort (NAKO).

The project investigates age-associated plasma protein signatures separately in men and women. The analysis pipeline includes data preprocessing and quality control, ordinal elastic-net feature selection, limma-based characterization of age associations, unsupervised clustering validation, protein network analysis, construction of the AgeProt score, and multiple-imputation sensitivity analyses.

## Repository structure

- `analysis/code/` – R scripts for preprocessing, main analyses, and sensitivity analyses
- `analysis/results/` – tables, figures, and analysis outputs

The main workflow is:

1. `01_preprocessing.R`
2. `02_analysis.R`
3. `03_mice_sensitivity.R`
4. `04_graphical_abstract.R`

## Data availability

Individual-level NAKO proteomic and participant data are not included in this repository due to data protection and access restrictions.

## Manuscript title
Sex-stratified plasma proteomic signatures of chronological aging in the Bremen study center of the German National Cohort
## Authors: 
Giovanni Corradini, Ronja Foraita, Stefan Loroch, Johannes Balkenhol, Kathrin Günther, Oliver Tüscher, Christian Herder, Tamara Schikowski, Helen Morrison, Krasimira Aleksandrova
