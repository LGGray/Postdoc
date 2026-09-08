#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_tile_floor_map
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

# spatial_tile_map.slurm with a depth floor, plus the distribution page that
# justifies the number. Same inputs, same resources, a separate output name -
# so the floored and unfloored figures sit side by side rather than one
# replacing the other.
#
#   sbatch slurm/spatial_tile_floor_map.slurm 64 9w,78w dup
#   sbatch slurm/spatial_tile_floor_map.slurm 64 9w,78w dup 50
#
# Read-only with respect to the tile run, but unlike spatial_tile_map.slurm it
# is meant for a FINISHED one: with a floor on, a below-floor tile and a
# not-yet-scored tile are both grey in the ratio panels. The script warns when
# it finds pending tiles rather than drawing the ambiguity silently.

module load slurm_setup

# Positional, not environment: --export=NONE stops `MIN_X_UNITS=20 sbatch` from
# reaching the job.
TILE_UM="${1:-64}"
SAMPLES="${2:-9w,78w}"
# Third arg is the directory suffix, as for spatial_tile_map.slurm: the SNP mask
# label for the original trees, or the filter mode ("dup") for the
# duplicate-inclusive one. no_Xist means "no suffix".
SNP_LABEL="${3:-no_Xist}"
# The floor, in the unit the tree was scored in. 20 is the ~20 independent
# observations panel_depth already names; at read level that is 5-8 molecules,
# so 50 is the honest read-level stand-in. The script's page 2 shows what
# either choice does to the pooled ratio.
MIN_X_UNITS="${4:-20}"
if [ "$SNP_LABEL" = "no_Xist" ]; then SUF=""; else SUF="_$SNP_LABEL"; fi

# Arg 3 doubles as the output-directory suffix, and for the duplicate-inclusive
# trees that suffix is a FILTER mode, not a SNP mask: spatial_sinto_tiles.slurm
# writes _dup and _raw while scoring against the no_Xist bed. Deriving the bed
# name from arg 3 would invent "SNPfile_..._dup.bed" and the provenance sidecar
# would carry no fingerprint, so map the filter suffixes back to the mask that
# was actually used. Arg 5 overrides, for a tile run scored against another.
case "$SNP_LABEL" in
  dup|raw|none) BED_DEFAULT=no_Xist ;;
  *)            BED_DEFAULT=$SNP_LABEL ;;
esac
SNP_BED_LABEL="${5:-$BED_DEFAULT}"

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
REPO=$BASE/Postdoc

if [ -z "$SCRATCH" ]; then
  echo "WARNING: SCRATCH is not set - only collected locus tables will be found" >&2
fi

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, NOT the repo-wide RNAseq env - see the note in
# spatial_tile_map.slurm. Do not "fix" this for consistency with the neighbours.
conda activate seurat_env

TILE_UM="$TILE_UM" SAMPLES="$SAMPLES" SNP_LABEL="$SNP_LABEL" \
  MIN_X_UNITS="$MIN_X_UNITS" SPATIAL_DIR="$REPO/spatial" \
  SNP_BED_LABEL="$SNP_BED_LABEL" \
  OUT_PDF="$BASE/adult_aged_spatial/ase/tile_ratio_map_${TILE_UM}um${SUF}_min${MIN_X_UNITS}.pdf" \
  Rscript "$REPO/spatial/tile_ratio_map_floor.R"

echo "Done. PDF, CSV and provenance in $BASE/adult_aged_spatial/ase/"
