#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J cellranger_arc
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
# ALTERNATIVE QUEUE. Swap back to these three lines - commenting out the serial
# two above - only once cm4 fairshare has recovered. serial_std takes no --qos;
# cm4_tiny does. See the long rationale in spatial_sinto_tiles.slurm: cm4_tiny
# was refusing to place jobs on FAIRSHARE, not on resources, so a smaller shape
# and more waiting both fail to help. cm4_tiny would be much the better home for
# this job on its merits (112 cpus / 256G against 16 / 96G here, so roughly a
# 6x wall-clock difference), so it is worth re-checking before a re-run:
#   sshare -U -u $USER    # has EffectvUsage decayed back toward NormShares?
#
##SBATCH --clusters=cm4
##SBATCH --partition=cm4_tiny
##SBATCH --qos=cm4_tiny
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=1-2

# 10x Multiome (GEX + ATAC) count for the two Eckstein/Engelhardt heart samples.
#
# Array 1 = 9w  (GEX 25L008363, ATAC 25L008352 isets 1-4)
# Array 2 = 78w (GEX 25L008364, ATAC 25L008356 isets 1-4)
#
# The ATAC "iset1..iset4" directories are the four i7 index oligos of ONE
# library, demultiplexed separately by the facility. They are listed as four
# fastqs rows of the same sample in libraries_*.csv -- do NOT treat them as
# four samples. Verified, not assumed: the top 1000 cell barcodes of iset1
# overlap isets 2/3/4 at 966/967/964, so all four sample the same nuclei.
#
# EXPECT THIS TO NEED MORE THAN ONE SUBMISSION. The GEX libraries are roughly
# 650-700M read pairs each (~5x the 10x recommendation, deliberate
# over-sequencing for allele-specific work) plus ~330M ATAC pairs, and
# serial_std caps us at 16 cpus and 24h. A rough estimate puts one sample at
# 20-35h, so budget two rounds. That is safe and cheap:
#
#   - cellranger-arc checkpoints per stage, so re-running the same --id in the
#     same directory resumes from the last completed stage.
#   - Resubmitting the WHOLE array is idempotent. The sample that already
#     finished detects a complete pipestance and exits in seconds; only the
#     unfinished one does work.
#
# So on a TIMEOUT, just `sbatch slurm/cellranger_arc_multiome.slurm` again.

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

# A killed pipestance leaves _lock behind, and cellranger-arc then refuses to
# restart rather than risk two writers. Since we expect to be killed by the 24h
# limit, clear it. This is only safe because one array task owns one --id: if
# you ever run the same ID twice CONCURRENTLY, this would let both proceed and
# corrupt the output.
if [ -f "$ID/_lock" ]; then
  echo "NOTE: clearing stale $ID/_lock left by a previous run, and resuming."
  rm -f "$ID/_lock"
fi

# --create-bam became mandatory in cellranger-arc 2.2.0 (2.0.2 defaulted it).
# It must be true here: every downstream allelic step reads BAMs, not the
# matrices -- sinto filterbarcodes (slurm/sinto.slurm), Allelome.PRO2
# (slurm/sc_Allelome.PRO2.slurm) and the pysam UMI route all need
# gex_possorted_bam.bam / atac_possorted_bam.bam. Budget the disk: at this
# depth expect roughly 100-120G of BAM per sample on top of the rest of outs/.
cellranger-arc count \
  --id="$ID" \
  --reference="$REF" \
  --libraries="$LIB" \
  --create-bam=true \
  --localcores=16 \
  --localmem=90
