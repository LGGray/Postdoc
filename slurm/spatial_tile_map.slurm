#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_tile_map
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00

# Draw the per-tile chrX allelic ratio onto the slide, for both samples, from
# whatever spatial_sinto_tiles.slurm has finished so far. Read-only: safe to run
# against a tile job that is still going.
#
#   sbatch ~/Postdoc/slurm/spatial_tile_map.slurm 64
#
# Not an array job - both samples are drawn in one pass so they land in one PDF,
# and the whole thing is 254 small text files plus tissue_positions.
#
# 16G/30min is generous. The cost is reading tissue_positions for two samples
# (~10M 2um bins each); the locus tables are trivial.

module load slurm_setup

# Positional, not environment: --export=NONE stops `TILE_UM=64 sbatch` from
# reaching the job. Must match the size the tile run used.
TILE_UM="${1:-64}"
SAMPLES="${2:-9w,78w}"

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
REPO=$BASE/Postdoc

# The script prefers the collected Allelome.PRO2_tiles_<um>um directory and
# falls back to $SCRATCH/spatial_tiles_<sample>_<um>um/allelome, which is where
# the tables sit until the tile job reaches its collect stage. Keep SCRATCH
# visible so the in-flight case works.
if [ -z "$SCRATCH" ]; then
  echo "WARNING: SCRATCH is not set - only collected locus tables will be found" >&2
fi

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, NOT the repo-wide RNAseq env. Most scripts in this directory use
# RNAseq and are misleading here; confirmed against seurat_env on 2026-08-27.
# Do not "fix" this for consistency with the neighbouring wrappers.
conda activate seurat_env

# arrow reads tissue_positions.parquet and png draws the H&E. Both are optional:
# the script falls back to a .tsv/.csv positions file and to no background
# image, so a missing package costs a panel, not the run.
Rscript -e 'for (p in c("arrow","png")) if (!requireNamespace(p, quietly=TRUE))
              message("note: R package \"", p, "\" not installed - see the script header")' \
  2>&1 | sed 's/^/  /'

TILE_UM="$TILE_UM" SAMPLES="$SAMPLES" \
  OUT_PDF="$BASE/adult_aged_spatial/ase/tile_ratio_map_${TILE_UM}um.pdf" \
  Rscript "$REPO/spatial/tile_ratio_map.R"

echo "Done. PDF and CSV in $BASE/adult_aged_spatial/ase/"
