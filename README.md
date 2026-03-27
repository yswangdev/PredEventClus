# Prediction of Event Times in Clinical Trials with Covariate-Driven Clustering

This repository contains the code and analysis scripts for the paper on predicting event times in clinical trials using covariate-driven clustering methods.

## Overview

This project implements and evaluates various clustering methods (PAM, Hierarchical, K-prototypes, Kamila) combined with survival analysis models (Exponential, Weibull, Gompertz) to predict event times in clinical trials. The approach accounts for patient heterogeneity through covariate-driven clustering to improve the predicting event times.


## Usage

### Simulation Studies

#### Scenario 1 (no true clustering; compare original vs clustering methods)
```r
source("Simulation/sim_by_no_clus_S1.R")
```

#### Scenario 2 (2 clusters, well separated)
```r
source("Simulation/sim_by_2clus_S2.R")
```

#### Scenario 3 (2 clusters, less separated)
```r
source("Simulation/sim_by_2clus_S3.R")
```

#### Scenario 4 (3 clusters, mid-level separated)
```r
source("Simulation/sim_by_3clus_S4.R")
```

### Real-World Data Analysis

```r
source("Real_World/real_data.R")
```

**Note:** This script requires pre-processed data files:
- `real_world_full_data.rds` - Full clinical trial dataset
- `real_world_mixed_data.rds` - Mixed-type covariates for clustering

### Visualization

Generate Kaplan-Meier curves:
```r
source("Visualization/KM_curve.R")
```

Generate MDS plots and other visualizations:
```r
source("Visualization/codes_for_plots.R")
```

## Output Files

### Simulation Results
- `predict_timepoints_no_clus_s1_n_[n].csv` / `predict_if_no_clus_s1_n_[n].csv` - Scenario 1 results
- `predict_timepoints_s2_n_[n].csv` / `predict_if_s2_n_[n].csv` - Scenario 2 results
- `predict_timepoints_s3_[scenario]_n_[n].csv` / `predict_if_s3_[scenario]_n_[n].csv` - Scenario 3 results
- `predict_timepoints_s4_n_[n].csv` / `predict_if_s4_n_[n].csv` - Scenario 4 results

### Real-World Results
- `[n]_[k]clus_BIC.csv` - Prediction metrics
- `[n]_[k]clus_cluster_labels.csv` - Cluster assignments


## Citation

If you use this code in your research, please cite the associated paper.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or issues, please contact the authors.
