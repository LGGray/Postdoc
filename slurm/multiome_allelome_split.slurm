#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J allelome_split
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
#SBATCH --mem=96G
#SBATCH --time=4:00:00

#   sbatch slurm/multiome_allelome_split.slurm 78w
#
# sinto split only, one sample, per nucleus. Separated from scoring so the
# scoring can then run as a job ARRAY: measured at 0.73 min/nucleus, 78w's 3125
# nuclei are 37.9h in one stream but 9.5h across four. The split itself cannot
# be parallelised the same way - sinto makes a single pass over the BAM
# whatever it is splitting - so it gets its own short job.
#
# Idempotent: skips if the BAM count already matches the number of distinct
# groups in the cells file.

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

ID="${1:?usage: sbatch multiome_allelome_split.slurm <9w|78w>}"
case "$ID" in 9w|78w) ;; *) echo "ID must be 9w or 78w" >&2; exit 1 ;; esac

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
WORK=$BASE/adult_aged_multiome
CELLS=$WORK/qc/sinto_${ID}_rna.txt
FLT=$WORK/allelome_pseudobulk/bam/${ID}_gex_dedup.bam
[ -f "$CELLS" ] || { echo "missing $CELLS" >&2; exit 1; }
[ -s "$FLT.bai" ] || { echo "missing filtered BAM $FLT (+.bai)" >&2; exit 1; }
if [ -z "${SCRATCH:-}" ]; then echo "SCRATCH is not set" >&2; exit 1; fi

BAMDIR=$SCRATCH/multiome_allelome/pernucleus/bam/$ID
mkdir -p "$BAMDIR"
nproc=$((SLURM_CPUS_PER_TASK - 2)); [ "$nproc" -lt 1 ] && nproc=1

# Count DISTINCT GROUPS, not lines. The cells file has one line per nucleus and
# a group column; for the per-nucleus files those coincide, but for a
# per-celltype file they do not, and comparing BAM count to line count made the
# run script re-split every time.
want=$(cut -f2 "$CELLS" | sort -u | wc -l)
have=$(ls "$BAMDIR"/*.bam 2>/dev/null | wc -l)
echo "$ID: $want groups wanted, $have BAMs present"
if [ "$have" -ge "$want" ]; then
  echo "split already complete, nothing to do"; exit 0
fi

conda activate sinto
echo "$(date +%T) splitting"
sinto filterbarcodes --bam "$FLT" --cells "$CELLS" --outdir "$BAMDIR" --nproc "$nproc"
echo "$(date +%T) done: $(ls "$BAMDIR"/*.bam 2>/dev/null | wc -l) BAMs"
