#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J arc_preflight
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00

# Pre-flight check on a 1M-read subsample, to be run BEFORE the full
# cellranger_arc_multiome.slurm job.
#
# Why this exists: the facility sequenced the GEX libraries with a symmetric
# 101/101 recipe instead of the 10x multiome 28/90 recipe. R1 is therefore
# 101bp rather than 28bp. The first 28bp are still barcode+UMI and
# cellranger-arc should read only those, but its fastq validator is stricter
# than cellranger's, so we confirm on a subsample rather than discovering a
# rejection 20 hours into a full-node job.
#
# HOW TO READ THE RESULT:
#   - Fails with a read-length / fastq-geometry error  -> geometry is the
#     problem; we must hard-trim GEX R1 to 28bp before running for real.
#   - Fails later, at cell calling / "insufficient cells" -> EXPECTED AND
#     GOOD. 1M reads cannot call cells. Getting that far means the geometry
#     was accepted and the full job is safe to launch.

module load slurm_setup

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
ARC_VER=2.2.0
export PATH="$BASE/cellranger-arc-$ARC_VER:$PATH"
REF="$BASE/refdata-cellranger-arc-GRCm39-2024-A"
[ -d "$REF" ] || { echo "missing reference $REF" >&2; exit 1; }

RAW=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/00_raw_data
ATAC_IN="$RAW/Project_1729_Eckstein_Engelhardt_multiomeHeart_snATAC/Sample_25L008352_iset1"
RNA_IN="$RAW/Project_1729_Eckstein_Engelhardt_multiomeHeart_snRNA/Sample_25L008363"

if [ -z "${SCRATCH:-}" ]; then echo "SCRATCH is not set - aborting" >&2; exit 1; fi
WORK="$SCRATCH/arc_preflight_${SLURM_JOB_ID}"
mkdir -p "$WORK/gex" "$WORK/atac" || exit 1
echo "Using scratch: $WORK"

NREADS=1000000
NLINES=$(( NREADS * 4 ))
P=He_9w_multiome_XBxC_XX

# head closes the pipe early, so zcat exits on SIGPIPE; that is fine here but
# would trip the pipefail we would otherwise want, hence no `set -o pipefail`.
# The 2>/dev/null above is only there for the SIGPIPE note. It must NOT be
# allowed to swallow a missing or unreadable input: that produces an empty
# subsample and cellranger-arc then runs for two hours on nothing. Test the
# file first.
for R in R1 R2; do
  f="$RNA_IN/${P}_S4_L001_${R}_001.fastq.gz"
  [ -r "$f" ] || { echo "ERROR: missing or unreadable GEX fastq: $f" >&2; exit 1; }
  zcat "$f" 2>/dev/null | head -n "$NLINES" \
    | gzip > "$WORK/gex/${P}_S4_L001_${R}_001.fastq.gz"
done
for R in R1 R2 R3; do
  f="$ATAC_IN/${P}_S13_L001_${R}_001.fastq.gz"
  [ -r "$f" ] || { echo "ERROR: missing or unreadable ATAC fastq: $f" >&2; exit 1; }
  zcat "$f" 2>/dev/null | head -n "$NLINES" \
    | gzip > "$WORK/atac/${P}_S13_L001_${R}_001.fastq.gz"
done

# And a short subsample is just as fatal as a missing one, so fail on it rather
# than printing it and carrying on.
echo "subsampled read counts (should each be $NREADS):"
short=0
for f in "$WORK"/gex/*.fastq.gz "$WORK"/atac/*.fastq.gz; do
  n=$(( $(zcat "$f" | wc -l) / 4 ))
  printf "  %-12s %s\n" "$n" "$(basename "$f")"
  [ "$n" -eq "$NREADS" ] || short=1
done
if [ "$short" -ne 0 ]; then
  echo "ERROR: at least one subsample is not $NREADS reads. The source fastq is" >&2
  echo "       shorter than that, or was truncated. Not running cellranger-arc." >&2
  exit 1
fi

{
  echo "fastqs,sample,library_type"
  echo "$WORK/gex,$P,Gene Expression"
  echo "$WORK/atac,$P,Chromatin Accessibility"
} > "$WORK/libraries_preflight.csv"
cat "$WORK/libraries_preflight.csv"

cd "$WORK" || exit 1
# Mandatory in 2.2.0. Kept true so this exercises the same code path as the
# real run, which needs BAMs for the allelic work.
cellranger-arc count \
  --id=preflight \
  --reference="$REF" \
  --libraries="$WORK/libraries_preflight.csv" \
  --create-bam=true \
  --localcores=16 \
  --localmem=60
