#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_marker_panels
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=4:00:00

# Marker panels for one or more Visium HD sections.
#
#   sbatch spatial_marker_panels.slurm            # 78w (default)
#   sbatch spatial_marker_panels.slurm 9w 78w     # both, one figure dir each
#
# Each sample writes to figures/<sample>_<bin>/, so re-running one age does not
# touch the other's panels. Walltime covers both sections; a single one takes
# well under an hour.

module load slurm_setup

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then SAMPLES=(78w); fi

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
BIN=${BIN:-square_016um}
cd $SPATIAL

# Fail before the queue time is spent rather than after loading the matrix.
for s in "${SAMPLES[@]}"; do
  d=$SPATIAL/$s/outs/binned_outputs/$BIN
  for f in "$d/filtered_feature_bc_matrix.h5" \
           "$d/analysis/umap/gene_expression_2_components/projection.csv" \
           "$d/analysis/clustering/gene_expression_graphclust/clusters.csv"; do
    if [ ! -s "$f" ]; then echo "Missing $f" >&2; exit 1; fi
  done
done

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, not RNAseq: the spatial R scripts are built against it, and this
# one reads tissue_positions.parquet, so it needs arrow the same way
# tile_ratio_map.R does. Do not "fix" this to RNAseq to match the neighbours.
conda activate seurat_env

Rscript -e '.libPaths(c("~/R/matrix-dev", .libPaths()));
            for (p in c("arrow", "patchwork", "ggrastr")) if (!requireNamespace(p, quietly=TRUE))
              message("note: R package \"", p, "\" not installed")' 2>&1 | sed 's/^/  /'

# SCTransform on ~60k 16um bins would be the memory peak, but the script runs
# LogNormalize by default and everything after it is plotting. Single-threaded
# on purpose - future::plan multisession duplicates the object per worker and
# this job is not slow enough to pay for that.
rc=0
for s in "${SAMPLES[@]}"; do
  echo "=== $s ($BIN) ==="
  SPATIAL_DIR="$SPATIAL" SAMPLE="$s" BIN="$BIN" \
    Rscript $BASE/Postdoc/spatial/spatial_marker_panels.R
  r=$?
  if [ $r -ne 0 ]; then echo "spatial_marker_panels.R failed for $s (exit $r)" >&2; rc=$r; fi
  echo
done

conda deactivate
exit $rc
