#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_v2_load_cells
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00

# Load the spaceranger cell-segmentation output into a Seurat object, one sample
# after another.
#
#   sbatch slurm/spatial_v2_load_cells.slurm            # 9w and 78w
#   sbatch slurm/spatial_v2_load_cells.slurm 9w         # one sample
#
# Positional, not environment: --export=NONE above means `MATRIX=... sbatch`
# does not reach the job. Override with --export=ALL,...:
#   sbatch --export=ALL,MATRIX=raw slurm/spatial_v2_load_cells.slurm 9w
#
# Resources: the filtered matrix is ~95k cells x 34k genes with 13.6M non-zeros,
# which is small; the peak is jsonlite holding the 50 MB cell geojson while the
# polygon areas are computed. 48G is generous for that and leaves room for the
# object and the figures. Each sample takes a few minutes. 4 cpus because
# nothing here threads - the figures are the slow part.

module load slurm_setup

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then SAMPLES=(9w 78w); fi

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
MATRIX=${MATRIX:-filtered}
cd "$SPATIAL" || exit 1

# Fail before the queue time is spent rather than after.
for s in "${SAMPLES[@]}"; do
  d=$SPATIAL/$s/outs/segmented_outputs
  for f in "$d/${MATRIX}_feature_cell_matrix.h5" "$d/cell_segmentations.geojson" \
           "$d/spatial/scalefactors_json.json"; do
    if [ ! -s "$f" ]; then echo "Missing $f" >&2; exit 1; fi
  done
  if [ ! -s "$d/nucleus_segmentations.geojson" ]; then
    echo "NOTE: no nucleus_segmentations.geojson for $s - nucleus columns will be NA" >&2
  fi
done

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, as for the spatial/ scripts: needs Seurat v5 (CreateFOV) and
# hdf5r to read the segmented h5.
conda activate seurat_env

Rscript -e '.libPaths(c("~/R/matrix-dev", .libPaths()));
            for (p in c("hdf5r", "jsonlite", "patchwork", "ggrastr", "data.table", "mgcv")) if (!requireNamespace(p, quietly=TRUE))
              message("note: R package \"", p, "\" not installed")' 2>&1 | sed 's/^/  /'

rc=0
for s in "${SAMPLES[@]}"; do
  echo "=== $s ($MATRIX) ==="
  SPATIAL_DIR="$SPATIAL" SAMPLE="$s" MATRIX="$MATRIX" \
    Rscript "$BASE/Postdoc/spatial_V2/01_load_segmented_cells.R"
  r=$?
  if [ $r -ne 0 ]; then echo "01_load_segmented_cells.R failed for $s (exit $r)" >&2; rc=$r; fi
  echo
done

conda deactivate
exit $rc
