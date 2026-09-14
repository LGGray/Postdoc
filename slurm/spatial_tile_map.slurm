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
#   sbatch ~/Postdoc/slurm/spatial_tile_map.slurm 64 9w,78w no_Xic500kb
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
# Which SNP mask the tile run used. Must match the label passed to
# spatial_sinto_tiles.slurm, since that is what named its output directory.
# no_Xist reads and writes the original paths; no_Xic500kb sits alongside them.
SNP_LABEL="${3:-no_Xist}"
if [ "$SNP_LABEL" = "no_Xist" ]; then SUF=""; else SUF="_$SNP_LABEL"; fi

# Arg 3 doubles as the output-directory suffix, and for the duplicate-inclusive
# trees that suffix is a FILTER mode, not a SNP mask: spatial_sinto_tiles.slurm
# writes _dup and _raw while scoring against the no_Xist bed. Deriving the bed
# name from arg 3 would invent "SNPfile_..._dup.bed" and the provenance sidecar
# would carry no fingerprint, so map the filter suffixes back to the mask that
# was actually used. Arg 4 overrides, for a tile run scored against another.
case "$SNP_LABEL" in
  dup|raw|none) BED_DEFAULT=no_Xist ;;
  *)            BED_DEFAULT=$SNP_LABEL ;;
esac
SNP_BED_LABEL="${4:-$BED_DEFAULT}"

# Arg 5: fold the ratio to max(A1,A2)/n (the ar_dom convention of
# OCM_heart/allelic_ratio) instead of A1/n. Positional like everything else,
# because --export=NONE means FOLD_RATIO=1 sbatch would not reach the job.
# It extends SUF so a folded run lands BESIDE the directional one rather than
# overwriting its PDF - OUT_PDF below is built from SUF and is the reason this
# cannot be a plain environment variable.
# The folded maps are a diagnostic, not a result: see the FOLD_RATIO note in
# spatial/tile_ratio_map.R for why folding is cheap on snRNA chrX (p~0.87) and
# expensive here (read-level pooled 0.60-0.72, where the fold adds a
# depth-dependent floor to a map whose depth varies across the section).
FOLD_RATIO="${5:-}"
case "$(echo "$FOLD_RATIO" | tr '[:upper:]' '[:lower:]')" in
  ""|0|false|no) FOLD_RATIO=""              ;;
  *)             FOLD_RATIO=1; SUF="${SUF}_folded" ;;
esac

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

TILE_UM="$TILE_UM" SAMPLES="$SAMPLES" SNP_LABEL="$SNP_LABEL" \
  SNP_BED_LABEL="$SNP_BED_LABEL" FOLD_RATIO="$FOLD_RATIO" \
  OUT_PDF="$BASE/adult_aged_spatial/ase/tile_ratio_map_${TILE_UM}um${SUF}.pdf" \
  Rscript "$REPO/spatial/tile_ratio_map.R"

echo "Done. PDF and CSV in $BASE/adult_aged_spatial/ase/"
