#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_mask_comparison
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
#SBATCH --mem=32G
#SBATCH --time=01:00:00

# Task 2 of spatial/NEXT_ANALYSIS.md: the per-allele rates under each SNP mask,
# side by side, so we know how much of the CAST rate difference survives masking
# the Xic rather than just the Xist gene span.
#
# Nothing is re-counted. The masked and unmasked columns both come out of the
# single counting pass in spatial_ase_sweep.slurm, which is why they are exactly
# comparable:
#
#   sbatch ~/Postdoc/slurm/spatial_ase_sweep.slurm 64 force
#   sbatch ~/Postdoc/slurm/spatial_mask_comparison.slurm 64
#
# If it reports that most of the effect was the Xic, the tile map has to be
# rebuilt on the no_Xic500kb SNP bed before it can be used - that is the
# expensive part, ~14h per sample:
#
#   sbatch ~/Postdoc/slurm/spatial_sinto_tiles.slurm 64 0 no_Xic500kb
#   sbatch ~/Postdoc/slurm/spatial_tile_map.slurm 64 9w,78w no_Xic500kb

module load slurm_setup

TILE_UM="${1:-64}"
SAMPLES="${2:-9w,78w}"
N_BOOT="${3:-10000}"

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
REPO=$BASE/Postdoc

missing=0
for s in ${SAMPLES//,/ }; do
  f=$SPATIAL/ase/$s/bin_allele_counts.tsv
  if [ ! -s "$f" ]; then echo "Missing: $f" >&2; missing=1; fi
done
if [ "$missing" -ne 0 ]; then
  echo "Run the counting pass first: sbatch $REPO/slurm/spatial_ase_sweep.slurm" >&2
  exit 1
fi

# Warn rather than fail: the script reports whatever sets are present and says
# so when the Xic-masked columns are absent.
for s in ${SAMPLES//,/ }; do
  if ! head -1 "$SPATIAL/ase/$s/bin_allele_counts.tsv" | grep -q 'noxic_ref'; then
    echo "NOTE: $s has no noxic_ref column, so there is no mask comparison to" >&2
    echo "      make. Build the interval beds and recount:" >&2
    echo "        cd $BASE/GRCm39 && OUT_DIR=. Rscript $REPO/OCM_heart/core_escape_SNPs.R" >&2
    echo "        sbatch $REPO/slurm/spatial_ase_sweep.slurm $TILE_UM force" >&2
  fi
done

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

BASE="$SPATIAL" TILE_UM="$TILE_UM" SAMPLES="$SAMPLES" N_BOOT="$N_BOOT" \
  Rscript "$REPO/spatial/ase_mask_comparison.R"

echo "Done. CSVs and PDF in $SPATIAL/ase/"
