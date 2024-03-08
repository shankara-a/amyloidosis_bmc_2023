# Subtypes of AL Amyloidosis
Github repository for manuscript titled:
    _Subtypes of AL Amyloidosis and Cumulative Incidence of Mortality and ESKD_

Work is completed with the Boston University Center for Amyloidosis.

Contact: Shankara Anand (shankara.k.anand@gmail.com)

## Repoistory Contents

All code to reproduce analysis + figures.

### Analysis Notebooks & Scripts
#### Data Processing
* `01a_process.ipynb`: process raw data, initial filtering analysis
* `01b_process_eskd.ipynb`: process ESKD treatment data

#### Data Imputation
* `02a_impute_mice.R`: run MICE imputation method for full dataset, 2004- dataset, and 2008 dataset
* `02b_impute_comparison.ipynb`: run KNN and median imputation methods for full dataset, 2004- dataset, and 2008 dataset; compute comparison metrics

#### Clustering
* `03a_cluster_ccp.R`: run consensusClusterPlus for full dataset, 2004- dataset, and 2008 dataset with missing values, MICE, KNN and median 
* `03b_cluster_comparison.ipynb`: analyzes all consensusClusterPlus runs, creates main figures
* `03c_cluster_ingroup.ipynb`: run statistical comparisons within select groups of interests (i.e. cardiac and renal stages); relevant for `Figure S4-5`
* `03c_FigS4-5.ipynb`: create alluvial plots for supplemental figures `S4-5`
* `03d_cluster_enrichment.R`: create balloonplots
* `03e_cluster_internal_validation.R`: run multiple clustering algorithms and compare result
* `03f_FigS3.ipynb`: create `Figure S3`

#### Survival Analyses
* `04a_adjsurv_mortality.ipynb`: notebook for all survival analysis (KM-curves, cox regression, adjusted survival curves) for all-cause mortality
* `04b_adjsurv_eskd.ipynb`: notebook for all survival analysis (KM-curves, cox regression, adjusted survival curves) for eskd

#### Prediction
* `05a_pred_subgroups.ipynb`: notebook for analyzing results of optimized prediction models
* `5b-d_pred_hyperopt*.py`: scripts for running bayesian hyperparameter optimization for all classification models for the prediction with baseline labs, with categorical clinical variables, and an abbreviated model with only top labs
* `5e*_pred_shap.py`: scripts for running shapley analysis and computed cross-validation metrics for all clusterings

#### Funcs
* `funcs/amyloid.py`: annotation information for dataset
* `funcs/plotting.py`: python plotting utils
* `funcs/utils.py`: analysis functions
* `funcs/rfuncs.R`: R functions for analyses

## Figures

### Main Figures
* `03b_cluster_comparison.ipynb`: 1A, 1C-E
* `03d_cluster_enrichment.R`: 2A-B
* `04a_adjsurv_mortality.ipynb`: 1B, 3A-B
* `04b_adjsurv_eskd.ipynb`: 4A-B
* `05a_pred_subgroups.ipynb`: 5A-D

### Supplemental Figures
* `01a_process.ipynb`: S1A, S1B, S1D
* `02b_impute_comparison.R`: S2A
* `03a_cluster_ccp.R`: S2B
* `03b_cluster_comparison.ipynb`: S1C, S1E, S2C-F
* `03c_FigS4-5.ipynb`: S4A, S5A
* `03d_cluster_enrichment.R`: S4B, S5B
* `03f_FigS3.ipynb`: S3
* `04a_adjsurv_mortality.ipynb`: S4C
* `04b_adjsurv_eskd.ipynb`: S5C

## Tables

### Main Tables
* `00_Table1-2.ipynb`: Table 1 + Table 2
* `04a_adjsurv_mortality.ipynb` + `04b_adjsurv_eskd.ipynb`: Table 3
* `05e1_pred_shap.py`: Table 4

### Supplemental Tables
* `02b_impute_comparison.ipynb`: Table S1-3
* `04a_adjsurv_mortality.ipynb`: Table S4
* `04b_adjsurv_eskd.ipynb`: Table S5
* `05e1_pred_shap.py`: Table S6