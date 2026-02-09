# SLE Classification Pipeline

Machine learning pipeline for classifying Systemic Lupus Erythematosus (SLE) from pseudobulk single-cell RNA-seq gene expression data using X chromosome and SLE-associated gene features.

## Overview

This pipeline implements a multi-model classification framework consisting of:

1. **Feature selection** (`boruta_enet.py`) — Boruta (Random Forest-based) and Elastic Net feature selection
2. **Model training** — Five classifiers trained with grid-searched hyperparameters:
   - Logistic Regression (`logit.regression.py`)
   - Random Forest (`random.forest.py`)
   - Support Vector Machine (`support.vector.machine.py`)
   - Gradient Boosting Machine (`gradient.boosting.machine.py`)
   - Multilayer Perceptron (`multilayer.perceptron.py`)
3. **Ensemble** (`voting.classifier.py`) — Soft voting ensemble of all five models
4. **Independent validation** (`predict_independent_SLE.py`) — Evaluation on an independent cohort
5. **Interpretability** (`SHAP.py`) — SHAP beeswarm plots for feature importance
6. **Visualisation** (`plot.all.curves.py`) — Overlaid PR curves for all models
7. **Statistical analysis** (R scripts) — Four focused analysis scripts:
   - `01_model_comparison.R` — Individual and ensemble model metric comparison, boxplots, and statistical tests
   - `02_feature_concordance.R` — Feature selection consistency across CV splits, heatmaps, cell composition (propeller), and escape gene enrichment
   - `03_independent_validation.R` — Independent cohort prediction metrics with confidence intervals
   - `04_chrX_feature_analysis.R` — chrX feature heatmaps, Jaccard similarity, DEG overlap, escapee analysis, and XIST/TSIX correlation

## Requirements

### Python
```
pip install -r requirements.txt
```

### R
The following R packages are required for the R analysis scripts and `rename_celltypes.R`:
- `dplyr`, `tidyr`, `tibble`, `stringr`, `reshape2`
- `ggplot2`, `ggpubr`, `ggsignif`, `gridExtra`, `ggrepel`
- `ComplexHeatmap`, `circlize`, `UpSetR`
- `Seurat`, `speckle`

## Setup

### Create directory structure

Run the setup script once from the project root to create all required directories:
```bash
./setup_directories.sh
```

This creates the full tree of ~290 directories across `pseudobulk/`, `pseudobulk_update/`, and `results_update/` needed by the pipeline.

## Usage

### 1. Feature Selection
```bash
python boruta_enet.py <input.RDS> <split_number>
```

### 2. Train Individual Models
```bash
python logit.regression.py <input.RDS> <split_number> <feature_type>
python random.forest.py <input.RDS> <split_number> <feature_type>
python support.vector.machine.py <input.RDS> <split_number> <feature_type>
python gradient.boosting.machine.py <input.RDS> <split_number> <feature_type>
python multilayer.perceptron.py <input.RDS> <split_number> <feature_type>
```

### 3. Ensemble Model
```bash
python voting.classifier.py <input.RDS> <split_number> <feature_type>
```

### 4. Independent Validation
```bash
python predict_independent_SLE.py <input.RDS> <split_number> <age_group>
```

### Arguments

| Argument | Description |
|----------|-------------|
| `input.RDS` | RDS file containing pseudobulk expression matrix with `class` and `individual` columns |
| `split_number` | Integer (1–10) specifying the cross-validation split |
| `feature_type` | Feature selection strategy: `intersection`, `combined`, `boruta`, or `enet` |
| `age_group` | For independent validation: `adult`, `child`, or `all` |

## Directory Structure

The pipeline expects the following directory layout:
```
pseudobulk_update/
  split_{1..10}/
    data.splits/       # Train/test CSV splits
    features/          # Selected features (Boruta, Elastic Net)
    scaler/            # StandardScaler pickle files
    {feature_type}/
      metrics/         # Classification metrics CSV files
      confusion/       # Confusion matrix PDFs
      AUROC/           # ROC curve PDFs
      PRC/             # PR curve PDFs
      ML.models/       # Trained model pickle files
      ensemble/        # Ensemble model outputs
results_update/        # Aggregated results and figures
```
