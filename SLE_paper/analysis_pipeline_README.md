# Analysis script

This directory contains a cleaned analysis script:

- `analysis_pipeline.R`

## What it does

- Loads training metrics across splits and methods
- Aggregates model performance metrics (including MCC)
- Runs statistical tests (`wilcox.test`, FDR correction)
- Produces figures and summary tables
- Loads feature selections and produces feature heatmaps
- Loads validation cohort metrics and generates validation figures

## Expected input files/directories

In the project root (or run directory), ensure these exist:

- `metadata_pbmc_female.control_managed.Rdata`
- `escapees.Rdata`
- `celltype.colours.RData`
- `X.immune.txt`
- `SLE_DisGeNet.tsv`
- `chrX_biomaRt.txt`
- `perez_update/` (split subfolders with `metrics`, `ensemble`, and `features`)
- `nehar_belaid_update/` (validation metrics)

## Run

From the project root:

```bash
Rscript SLE_paper/analysis_pipeline.R
```

## Outputs

Outputs are written to:

- `analysis_outputs/`

Key files include:

- `models_metrics_df.RData`
- `chrX_vs_SLE_wilcox_results.csv`
- `compare_models_mcc_combined.pdf`
- `compare_avg_mcc_combined_X_SLE.pdf`
- `all_celltypes_top_features_heatmap.pdf`
- `all_celltypes_all_features_heatmap.pdf`
- `validation_mcc.pdf`
