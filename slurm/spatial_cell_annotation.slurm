#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_cell_annotation
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00

# Cell-type annotation of the Visium HD sections at 8um bins, one sample after
# another, then the cross-sample composition figure once both exist.
#
#   sbatch spatial_cell_annotation.slurm            # 9w and 78w
#   sbatch spatial_cell_annotation.slurm 78w        # one sample
#
# Positional, not environment: --export=NONE above means `BIN=... sbatch` does
# not reach the job. Override the bin or the transfer with --export=ALL,...:
#   sbatch --export=ALL,BIN=square_016um,DO_TRANSFER=0 spatial_cell_annotation.slurm 9w
#
# Memory: the 8um matrix is ~250k bins; label transfer sketches to 50k before
# anchoring, and the OCM reference object is loaded alongside. 96G is the
# serial_std ceiling used elsewhere in this project and has headroom. Each
# sample is well under two hours with transfer, minutes without it.

module load slurm_setup

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then SAMPLES=(9w 78w); fi

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
BIN=${BIN:-square_008um}
cd "$SPATIAL" || exit 1

# Fail before the queue time is spent rather than after loading the matrix.
for s in "${SAMPLES[@]}"; do
  d=$SPATIAL/$s/outs/binned_outputs/$BIN
  for f in "$d/filtered_feature_bc_matrix.h5"; do
    if [ ! -s "$f" ]; then echo "Missing $f" >&2; exit 1; fi
  done
  if [ ! -d "$SPATIAL/$s/outs/binned_outputs/square_002um/spatial" ]; then
    echo "NOTE: no square_002um for $s - domain_2um.tsv will not be written" >&2
  fi
done
REF=$BASE/OCM/heart_seurat_object_SCT.rds
if [ ! -s "$REF" ]; then
  echo "NOTE: no snRNA reference at $REF - label transfer will be skipped" >&2
fi

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, as for spatial_marker_panels.R: needs Seurat v5 (SketchData,
# ProjectData) and arrow for tissue_positions.parquet.
conda activate seurat_env

Rscript -e '.libPaths(c("~/R/matrix-dev", .libPaths()));
            for (p in c("arrow", "patchwork", "ggrastr", "data.table")) if (!requireNamespace(p, quietly=TRUE))
              message("note: R package \"", p, "\" not installed")' 2>&1 | sed 's/^/  /'

rc=0
for s in "${SAMPLES[@]}"; do
  echo "=== $s ($BIN) ==="
  SPATIAL_DIR="$SPATIAL" SAMPLE="$s" BIN="$BIN" REF_RDS="$REF" \
    Rscript "$BASE/Postdoc/spatial/spatial_cell_annotation.R"
  r=$?
  if [ $r -ne 0 ]; then echo "spatial_cell_annotation.R failed for $s (exit $r)" >&2; rc=$r; fi
  echo
done

conda deactivate
exit $rc
