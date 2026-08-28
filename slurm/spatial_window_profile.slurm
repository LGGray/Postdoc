#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_window_profile
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

# Task 1 of spatial/NEXT_ANALYSIS.md: locate the excess chrX UMIs on the
# chromosome. Decides whether the +45% CAST per-allele rate is a locus artefact
# or a genuine Xi-wide effect, and therefore whether tasks 5-8 are worth running.
#
#   sbatch ~/Postdoc/slurm/spatial_ase_sweep.slurm 64 force   # writes the input
#   sbatch ~/Postdoc/slurm/spatial_window_profile.slurm 64    # this job
#
# Not an array job: both samples are analysed together, because the comparison
# between them is half the output.
#
# 32G/1h is generous. The cost is the two per-2um-bin count tables (a few
# million rows each); the window tables are sparse and small.

module load slurm_setup

# Positional, not environment: --export=NONE stops `TILE_UM=64 sbatch` from
# reaching the job. MUST match the --window-tile-um the counting pass used, or
# the tile grids do not line up - the script checks and refuses rather than
# plotting an empty panel.
TILE_UM="${1:-64}"
SAMPLES="${2:-9w,78w}"

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
REPO=$BASE/Postdoc

missing=0
for s in ${SAMPLES//,/ }; do
  w=$SPATIAL/ase/$s/chrX_windows_${TILE_UM}um.tsv.gz
  c=$SPATIAL/ase/$s/bin_allele_counts.tsv
  for f in "$w" "$c"; do
    if [ ! -s "$f" ]; then echo "Missing: $f" >&2; missing=1; fi
  done
done
if [ "$missing" -ne 0 ]; then
  echo "" >&2
  echo "The window table comes from the counting pass. Re-run it with:" >&2
  echo "  sbatch $REPO/slurm/spatial_ase_sweep.slurm ${TILE_UM} force ${TILE_UM}" >&2
  exit 1
fi

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

BASE="$SPATIAL" TILE_UM="$TILE_UM" SAMPLES="$SAMPLES" \
  Rscript "$REPO/spatial/ase_chrX_window_profile.R"

echo "Done. PDF and CSVs in $SPATIAL/ase/"
