#!/usr/bin/env bash
# setup_directories.sh
#
# Creates the full directory structure required by the SLE classification
# pipeline. Run once from the project root before executing the pipeline.
#
# Usage:
#   chmod +x setup_directories.sh
#   ./setup_directories.sh

set -euo pipefail

SPLITS=$(seq 1 10)
FEATURE_TYPES="boruta enet intersection combined"
MODEL_SUBDIRS="metrics confusion AUROC PRC ML.models ensemble"

echo "Creating SLE pipeline directory structure..."

# --- Input split indices (one-time, user provides these) ---
for s in $SPLITS; do
    mkdir -p "pseudobulk/split_${s}"
done

# --- Main working tree ---
for s in $SPLITS; do
    # Shared directories per split
    mkdir -p "pseudobulk_update/split_${s}/data.splits"
    mkdir -p "pseudobulk_update/split_${s}/features"
    mkdir -p "pseudobulk_update/split_${s}/scaler"
    mkdir -p "pseudobulk_update/split_${s}/feature.select.model"

    # Per feature-type directories
    for ft in $FEATURE_TYPES; do
        for sub in $MODEL_SUBDIRS; do
            mkdir -p "pseudobulk_update/split_${s}/${ft}/${sub}"
        done
    done
done

# --- Results / output directories ---
mkdir -p results_update/SHAP
mkdir -p results_update/feature_heatmap
mkdir -p results_update/jaccard_heatmaps

echo "Done. Directory structure created:"
echo ""
# Print a compact summary
echo "  pseudobulk/split_{1..10}/"
echo "  pseudobulk_update/split_{1..10}/"
echo "    ├── data.splits/"
echo "    ├── features/"
echo "    ├── scaler/"
echo "    ├── feature.select.model/"
echo "    └── {boruta,enet,intersection,combined}/"
echo "        ├── metrics/"
echo "        ├── confusion/"
echo "        ├── AUROC/"
echo "        ├── PRC/"
echo "        ├── ML.models/"
echo "        └── ensemble/"
echo "  results_update/"
echo "    ├── SHAP/"
echo "    ├── feature_heatmap/"
echo "    └── jaccard_heatmaps/"
echo ""
echo "Total directories created: $(find pseudobulk pseudobulk_update results_update -type d | wc -l | tr -d ' ')"
