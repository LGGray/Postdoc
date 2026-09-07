#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J cellranger_arc
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --array=1-2

# 10x Multiome (GEX + ATAC) count for the two Eckstein/Engelhardt heart samples.
#
# Array 1 = 9w  (GEX 25L008363, ATAC 25L008352 isets 1-4)
# Array 2 = 78w (GEX 25L008364, ATAC 25L008356 isets 1-4)
#
# The ATAC "iset1..iset4" directories are the four index oligos of ONE sample
# index set, demultiplexed separately by the facility. They are a single
# library and are listed as four fastqs rows of the same sample in
# libraries_*.csv -- do NOT treat them as four samples.
#
# Depth note: the GEX libraries are ~650-700M read pairs each (~5x the 10x
# recommendation), deliberately over-sequenced for allele-specific work. That
# is what drives the full-node request below. If 24h proves short, just
# resubmit: cellranger-arc checkpoints per stage and re-running the same --id
# in the same directory resumes rather than restarting.

module load slurm_setup

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2

# Edit ARC_VER if you install a different version.
ARC_VER=2.2.0
export PATH="$BASE/cellranger-arc-$ARC_VER:$PATH"
REF="$BASE/refdata-cellranger-arc-GRCm39-2024-A"

WORK="$BASE/adult_aged_multiome"
mkdir -p "$WORK"
cd "$WORK" || { echo "missing $WORK" >&2; exit 1; }

case "${SLURM_ARRAY_TASK_ID:-0}" in
  1) ID=9w  ;;
  2) ID=78w ;;
  *) echo "unexpected SLURM_ARRAY_TASK_ID '${SLURM_ARRAY_TASK_ID:-}'" >&2; exit 1 ;;
esac

LIB="$BASE/Postdoc/multiome/libraries_${ID}.csv"
[ -f "$LIB" ] || { echo "missing $LIB" >&2; exit 1; }
[ -d "$REF" ] || { echo "missing reference $REF" >&2; exit 1; }

cellranger-arc count \
  --id="$ID" \
  --reference="$REF" \
  --libraries="$LIB" \
  --localcores=110 \
  --localmem=250
