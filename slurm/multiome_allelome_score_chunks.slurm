#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J allelome_score
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-4

#   sbatch slurm/multiome_allelome_score_chunks.slurm 78w
#
# Score per-nucleus BAMs in parallel slices. Requires the split to exist -
# run multiome_allelome_split.slurm first.
#
# WHY AN ARRAY. Measured on 849 nuclei of 9w: 0.73 min/nucleus at 14-way
# parallel within one job. 78w's 3125 nuclei are therefore 37.9h in a single
# stream, which is three nights against a 24h limit, but 9.5h across four
# tasks. The serial cluster caps a user at 96 concurrent cores, so 4 tasks x 16
# leaves room alongside another 16-cpu job; do not raise the array width past 5
# without accounting for whatever else is running.
#
# Slices are taken by round-robin over the sorted BAM list, so every task gets
# a mix rather than one task inheriting all the large nuclei.
#
# Resume is per nucleus, so a task that is killed loses only what it had in
# flight, and a stale output directory left by a killed run is REMOVED before
# retrying - Allelome.PRO2 derives its output directory name from the date and
# the BAM name, so a half-written directory from the same day would otherwise
# collide with the retry.

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq
module load parallel
shopt -s nullglob

ID="${1:?usage: sbatch multiome_allelome_score_chunks.slurm <9w|78w>}"
case "$ID" in 9w|78w) ;; *) echo "ID must be 9w or 78w" >&2; exit 1 ;; esac

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
GRCM39=$BASE/GRCm39
PIPELINE=$BASE/Allelome.LINK/00_src/Allelome.PRO2.sh
ANNOTATION=$GRCM39/chr_annotation_mm39.bed
SNPS=$GRCM39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed
if [ -z "${SCRATCH:-}" ]; then echo "SCRATCH is not set" >&2; exit 1; fi
SCR=$SCRATCH/multiome_allelome
BAMDIR=$SCR/pernucleus/bam/$ID
OUTDIR=$SCR/pernucleus/out/$ID
mkdir -p "$OUTDIR"

N=${SLURM_ARRAY_TASK_COUNT:-1}
I=${SLURM_ARRAY_TASK_ID:-1}
nproc=$((SLURM_CPUS_PER_TASK - 2)); [ "$nproc" -lt 1 ] && nproc=1

n_bam=$(ls "$BAMDIR"/*.bam 2>/dev/null | wc -l)
[ "$n_bam" -gt 0 ] || { echo "no BAMs in $BAMDIR - run multiome_allelome_split.slurm $ID first" >&2; exit 1; }
echo "task $I/$N over $n_bam BAMs in $BAMDIR, -j $nproc"

MINE=$SCR/pernucleus/slice_${ID}_${I}of${N}.txt
ls "$BAMDIR"/*.bam | sort | awk -v n="$N" -v i="$I" 'NR % n == i % n' > "$MINE"
echo "slice size: $(wc -l < "$MINE")"

TASKS=$SCR/pernucleus/tasks_${ID}_${I}of${N}.txt; : > "$TASKS"
done_n=0; todo_n=0; stale_n=0
while IFS= read -r b; do
  base=$(basename "$b" .bam)
  if compgen -G "$OUTDIR/*_${base}_*/locus_table.txt" > /dev/null; then
    done_n=$((done_n+1)); continue
  fi
  # remove a half-written directory from a killed run before retrying
  for d in "$OUTDIR"/*_"${base}"_*/; do
    [ -d "$d" ] && { rm -rf "$d"; stale_n=$((stale_n+1)); }
  done
  echo "bash $PIPELINE -i '$b' -a $ANNOTATION -s $SNPS -r 1 -t 1 -o '$OUTDIR'" >> "$TASKS"
  todo_n=$((todo_n+1))
done < "$MINE"
echo "already scored: $done_n   to do: $todo_n   stale dirs cleared: $stale_n"

[ "$todo_n" -eq 0 ] && { echo "nothing to do"; exit 0; }
echo "$(date +%T) scoring"
parallel -j "$nproc" :::: "$TASKS"
echo "$(date +%T) task $I done"
