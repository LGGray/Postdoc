#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J sc_Allelome.PRO2.dedup
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --cpus-per-task=110
#SBATCH --mem=150G
#SBATCH --ntasks=1
#SBATCH --qos=cm4_tiny
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --array=1-4

# Allelome.PRO2 over the per-cell snRNA-seq BAMs, genome-wide, in either of two
# modes. FILTER=1 (default) scores UMI-collapsed, uniquely-mapped reads; FILTER=0
# scores every read, the way scAllelome.slurm did.
#
#   sbatch scAllelome_dedup.slurm              # -> Allelome.PRO2_dedup/
#   FILTER=0 sbatch scAllelome_dedup.slurm     # -> Allelome.PRO2/
#
# Two modes in one script rather than two scripts, so the only difference
# between the trees is the filter. Everything else - annotation, SNP file,
# thresholds, resume logic - is shared by construction and cannot drift.
#
# Why the trees have to coexist rather than one replacing the other:
# `total_reads` in locus_table.txt changes units between them, from reads to
# (approximately) molecules. Every downstream threshold calibrated on the
# read-level tree - MIN_TOTAL_READS in 02_whole_chrX.R, the sweep grids in 06/07,
# the exact-count stratification in 08 - is not transferable. 09 compares them
# cell by cell before anything downstream is re-derived.
#
# Measured on 50 9w cells (1.14M reads): 39.4% flagged duplicate, 7.4% not
# uniquely mapped (MAPQ 0/1/3 in STAR's encoding, of which MAPQ 0 - five or more
# loci - is the largest class). Expect total_reads to roughly halve.
#
#   -F 0x400  drop PCR duplicates. Cellranger marks these per (CB, UB, gene), so
#             what survives is ~one read per molecule. Every read of one UMI
#             carries the same allele, so duplicates add weight to a cell's ratio
#             without adding evidence about which allele was there.
#   -q 255    STAR-unique only. This is the post-hoc equivalent of
#             --outFilterMultimapNmax 1 in the bulk mapping scripts. It matters
#             more here than for bulk expression: a multimapper overlapping a SNP
#             is not merely misattributed, it is assigned an allele outright.
#
# Both filters match ase_bin_allele_counts.py and spatial_sinto_tiles.slurm, so
# the snRNA-seq and spatial numbers become comparable.

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM

nproc=$((SLURM_CPUS_PER_TASK - 2))
if [ "$nproc" -lt 1 ]; then nproc=1; fi

module load parallel
shopt -s nullglob

PIPELINE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh
SNPS=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed
GRCM39=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39
MIN_READ=1
TOTAL_READ=1
FILTER_MAPQ=255

# Positional, not environment: --export=NONE above stops `FILTER=0 sbatch ...`
# from reaching the job. The environment is kept as a fallback for interactive
# use only.
FILTER="${1:-${FILTER:-1}}"
case "$FILTER" in
  1) TREE=Allelome.PRO2_dedup ;;
  0) TREE=Allelome.PRO2 ;;
  *) echo "FILTER must be 0 or 1, got '$FILTER'" >&2; exit 1 ;;
esac
echo "FILTER=$FILTER -> $TREE"

for f in "$PIPELINE" "$SNPS"; do
  [ -f "$f" ] || { echo "Missing file: $f" >&2; exit 1; }
done

run_sample () {
  local SAMPLE=$1 ANNOTATION=$2
  local SRC=sinto/$SAMPLE
  local FLT=sinto_flt/$SAMPLE
  local OUT=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/$TREE/$SAMPLE
  local ANNOT_BASE SCORE_DIR
  ANNOT_BASE=$(basename "$ANNOTATION")

  if [ ! -d "$SRC" ]; then
    echo "SKIP $SAMPLE: $SRC does not exist" >&2; return 0
  fi
  [ -f "$ANNOTATION" ] || { echo "Missing annotation: $ANNOTATION" >&2; return 1; }
  mkdir -p "$OUT"

  # Refuse to add to a tree that already holds results from a DIFFERENT
  # annotation. Allelome.PRO2 puts the bed name in each run directory, so two
  # annotations can sit side by side in one sample directory and the loader in
  # 00_functions.R will read both - giving a cell one row per gene from the old
  # bed and one row per chromosome from the new one. That is precisely how 9w
  # ended up being scored a gene at a time. Fail here instead.
  local foreign
  foreign=$(ls -d "$OUT"/*_*.bed_*/ 2>/dev/null | grep -v "_${ANNOT_BASE}_" | head -3)
  if [ -n "$foreign" ]; then
    echo "ERROR: $OUT already holds runs from another annotation:" >&2
    echo "$foreign" >&2
    echo "       Move them aside before scoring with $ANNOT_BASE, e.g." >&2
    echo "       mkdir -p ${OUT}_pergene && mv $OUT/*_annotation_us_*.bed_*/ ${OUT}_pergene/" >&2
    return 1
  fi

  ##### ---- filter duplicates and multimappers ---- #####
  if [ "$FILTER" = "1" ]; then
    mkdir -p "$FLT"
    # Resume on the index, not the BAM: a job killed mid-`samtools view` leaves
    # a truncated but non-empty BAM, and the `&&` means it never got an index.
    # The index is the only marker that says the pair of steps finished.
    local list=$OUT/filter-list.txt
    : > "$list"
    for bam in "$SRC"/*.bam; do
      local c; c=$(basename "$bam")
      [ -s "$FLT/$c.bai" ] && continue
      echo "samtools view -b -F 0x400 -q $FILTER_MAPQ -o '$FLT/$c' '$bam' && samtools index '$FLT/$c'"
    done >> "$list"

    local n_flt; n_flt=$(wc -l < "$list")
    if [ "$n_flt" -gt 0 ]; then
      echo "[$(date +%T)] $SAMPLE: filtering $n_flt cell BAMs"
      parallel -j "$nproc" :::: "$list"
    else
      echo "[$(date +%T)] $SAMPLE: filtered cell BAMs already present"
    fi

    local ex exb
    ex=$(ls "$SRC"/*.bam | head -1); exb=$(basename "$ex")
    echo "  example $exb: $(samtools view -c "$ex") -> $(samtools view -c "$FLT/$exb") reads"
    SCORE_DIR=$FLT
  else
    SCORE_DIR=$SRC
    echo "[$(date +%T)] $SAMPLE: FILTER=0, scoring raw cell BAMs"
    echo "  total_reads in this tree counts reads, not molecules; ~39% are PCR"
    echo "  duplicates and carry no independent allele evidence."
  fi

  ##### ---- Allelome.PRO2, one invocation per cell ---- #####
  # Allelome.PRO2 writes a dated job directory per run,
  #   <out>/2026_08_27_<bam basename>_<annotation basename>_<run no.>/
  # and the trailing counter increments on every re-run. Without a resume check
  # a resubmitted job silently accumulates _2, _3, ... copies of finished work.
  local tasks=$OUT/task-list.txt
  : > "$tasks"
  local n_todo=0 n_done=0
  for bam in "$SCORE_DIR"/*.bam; do
    local cell; cell=$(basename "$bam" .bam)
    if compgen -G "$OUT/*_${cell}_${ANNOT_BASE}_*/locus_table.txt" > /dev/null; then
      n_done=$((n_done + 1)); continue
    fi
    echo "$PIPELINE -i $SCORE_DIR/$cell.bam -a $ANNOTATION -s $SNPS -r $MIN_READ -t $TOTAL_READ -o $OUT"
    n_todo=$((n_todo + 1))
  done >> "$tasks"

  echo "[$(date +%T)] $SAMPLE: $n_todo to run, $n_done already done"
  [ "$n_todo" -eq 0 ] && return 0

  local start=$SECONDS
  parallel -j "$nproc" --joblog "$OUT/parallel.log" :::: "$tasks"
  echo "[$(date +%T)] $SAMPLE: $n_todo cells in $((SECONDS - start))s at -j $nproc"

  local failed; failed=$(awk 'NR > 1 && $7 != 0' "$OUT/parallel.log" | wc -l)
  echo "  Allelome.PRO2 failures this pass: $failed / $n_todo"
  echo "  locus tables now present: $(find "$OUT" -name locus_table.txt | wc -l)"
}

# Genome-wide for all four, including 9w. scAllelome.slurm gave 9w a chrX-only
# per-gene bed (annotation_us_mm39_chrX.bed, 3227 intervals) and the other three
# a whole-chromosome one (chr_annotation_mm39.bed, 21 intervals). Two things
# went wrong as a result: 9w had no autosomal rows at all, hence no technical
# control, and 3227 chrX rows per cell meant the whole-chromosome match() in 02
# and 06 silently handed each cell one arbitrary gene. See
# assert_one_row_per_cell() in 00_functions.R.
#
# Running FILTER=0 here rebuilds the read-level tree for 9w on the same
# annotation as everything else, so all four samples are treated identically in
# both trees and 09 can pair every sample. The existing per-gene 9w results have
# to be moved aside first - the guard in run_sample refuses to mix annotations
# in one directory and prints the command.
#
# One array task per sample rather than four in sequence. Each sample gets its
# own 24h wall clock instead of sharing one, and a task that times out is
# resubmitted on its own:
#
#   sbatch scAllelome_dedup.slurm              # dedup tree, all four
#   sbatch scAllelome_dedup.slurm 0            # read-level tree, all four
#   sbatch --array=1 scAllelome_dedup.slurm 0  # read-level, 9w only
#
# Resume is per sample and keyed on the locus tables, so a resubmitted task
# picks up where it stopped rather than starting over, and FILTER=0 over
# 78w/Sham/TAC is a no-op that finds every cell already done. Add %2 to the
# array spec if four whole nodes queue badly.
case "$SLURM_ARRAY_TASK_ID" in
  1) SAMPLE=9w   ;;
  2) SAMPLE=78w  ;;
  3) SAMPLE=Sham ;;
  4) SAMPLE=TAC  ;;
  *) echo "Unexpected SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

# Propagate the failure. Without this the array task exits 0 on a missing
# annotation and the mail says END rather than FAIL.
run_sample "$SAMPLE" "$GRCM39/chr_annotation_mm39.bed" || exit 1

echo "[$(date +%T)] $SAMPLE done"
