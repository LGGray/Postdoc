#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_allelome_atac
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
#SBATCH --time=24:00:00

# Allelome.PRO2 over the ATAC modality of the same nuclei, mirroring
# multiome_allelome_run.slurm stage for stage so the two are directly
# comparable. Step 7 of multiome/ANALYSIS_PLAN.md.
#
# WHAT THIS MEASURES, AND WHY IT IS NOT THE SAME NUMBER AS THE RNA RUN.
# CAST is the inactive X in every nucleus (Xist is deleted on B6), so on the
# GEX BAM the CAST fraction on chrX is escape TRANSCRIPTION, ~12.7%. Here it is
# ACCESSIBILITY of the inactive X - a different quantity, and expected to be
# lower than the RNA escape wherever a locus is open but not transcribed. The
# two are only comparable as a concordance, never as the same number.
#
# CHROMOSOME-LEVEL ANNOTATION, NOT PEAKS, DELIBERATELY. chr_annotation_mm39.bed
# is the identical annotation the GEX run used, so this pass answers the gate
# question (is there Xi accessibility at all, and what is the autosomal mapping
# bias for 50bp BWA reads) against a directly comparable denominator and with
# no new assets. Peaks are the follow-up and need a decision this script does
# not make: 9w called 52,231 peaks and 78w 71,694, so a per-peak age contrast
# needs a union peak set rather than each sample's own.
#
# ---------------------------------------------------------------------------
# THE THREE THINGS THAT DIFFER FROM THE GEX RUN
# ---------------------------------------------------------------------------
#
#  1. -q 30, NOT -q 255. GEX is STAR, where unique = 255. ATAC is BWA, which
#     caps at 60, so -q 255 returns ZERO reads without erroring. The filter is
#     not cosmetic: the first chrX read sampled from this BAM is MAPQ 5 with an
#     XA tag listing four alternative placements, and mpileup applies no MAPQ
#     filter of its own (Allelome.PRO2.sh:238), so without -q 30 that read is
#     assigned an allele outright.
#
#  2. -F 0x400 BUYS MUCH LESS HERE. On GEX, 98.58%/97.79% duplicates meant a
#     ~70x smaller BAM. ATAC is 77.74%/66.93% duplicate, so expect ~3-4x. It is
#     still a no-op for CORRECTNESS either way - mpileup's default exclusion
#     mask already drops DUP - so this is purely to make sinto and the pileup
#     tractable. Note ATAC has no UMIs: deduplication is positional only, which
#     is cellranger-arc's call, not something this script can improve on.
#
#  3. THE CELLS FILES ARE THE `_rna` ONES. This looks wrong and is not. See
#     below - it is the single thing most likely to be "fixed" into silence.
#
# ---------------------------------------------------------------------------
# BARCODE SPACE - MEASURED, AND THE OPPOSITE OF WHAT 03 ASSUMED
# ---------------------------------------------------------------------------
# cellranger-arc keys the two modalities differently and per_barcode_metrics.csv
# carries all three columns (barcode, gex_barcode, atac_barcode). The ATAC BAM
# was assumed to be atac_barcode-keyed, and 03_signac_joint.R exported
# sinto_<id>_atac_bycelltype.txt on that basis. Sampling the BAM says otherwise:
#
#   CR:Z:TTGCTCCGTTAGCTCG      raw ATAC barcode
#   CB:Z:GACTCACCATGACTAT-1    corrected - and this is the CANONICAL barcode
#
# and per_barcode_metrics.csv puts exactly that pair on one row. Over 25 distinct
# CB tags sampled from three loci: 25/25 match gex_barcode, 0/25 match
# atac_barcode. cellranger-arc has already translated CB; CR holds the raw one.
#
# So the ATAC BAM takes the SAME cells files as the GEX BAM, and the
# sinto_<id>_atac_bycelltype.txt files are the ones that would produce EMPTY
# BAMs - silently, since sinto writes an empty file rather than failing on a
# barcode that matches nothing. That is why stage 0 below is a hard preflight
# rather than a comment: this failure mode has no error message.
#
# ---------------------------------------------------------------------------
# STAGES, cheapest first, EVERY ONE RESUMABLE
# ---------------------------------------------------------------------------
#   0  preflight: prove CB overlaps the cells file        <- seconds
#   A  filter both ATAC BAMs once   (-F 0x400 -q 30)      <- ~1h
#   B  PER CELL TYPE - 7 groups x 2 samples, minutes      <- guaranteed to land
#   C  PER NUCLEUS   - 4547 BAMs, hours to days           <- the actual question
#
# Each unit is skipped when its output exists, so a TIMEOUT plus a resubmit
# continues rather than restarting. Expect to resubmit: the post-filter ATAC BAM
# is roughly 15x the size of the GEX one, and stage C on GEX already took hours.
# If it times out twice, switch to serial_long (168h, --qos=cm4_serial_long,
# which caps the user at 100G across ALL running jobs, so run nothing else).
#
#   sbatch slurm/multiome_allelome_atac.slurm
# ---------------------------------------------------------------------------

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
GRCM39=$BASE/GRCm39
WORK=$BASE/adult_aged_multiome
PIPELINE=$BASE/Allelome.LINK/00_src/Allelome.PRO2.sh
ANNOTATION=$GRCM39/chr_annotation_mm39.bed
SNPS=$GRCM39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed
FIGS=$WORK/figures_signac        # per-celltype cells files from 03
QC=$WORK/qc                      # per-nucleus cells files from 01
SAMPLES="9w 78w"

for f in "$PIPELINE" "$ANNOTATION" "$SNPS"; do
  [ -f "$f" ] || { echo "Missing: $f" >&2; exit 1; }
done
if [ -z "${SCRATCH:-}" ]; then echo "SCRATCH is not set - aborting" >&2; exit 1; fi

SCR=$SCRATCH/multiome_allelome_atac
mkdir -p "$SCR"
nproc=$((SLURM_CPUS_PER_TASK - 2)); [ "$nproc" -lt 1 ] && nproc=1
module load parallel
shopt -s nullglob
echo "scratch: $SCR   parallel jobs: $nproc"

conda activate RNAseq

##### ============ stage 0: preflight, barcode space ============ #####
# The one check worth its own stage. A cells file keyed in the wrong barcode
# space costs a full run and reports success, so refuse to start without
# evidence that CB and the cells file are in the same space.
echo
echo "[0] $(date +%T) preflight: barcode space"
for ID in $SAMPLES; do
  SRC=$WORK/$ID/outs/atac_possorted_bam.bam
  CELLS=$QC/sinto_${ID}_rna.txt
  [ -f "$SRC" ]   || { echo "[0] missing $SRC" >&2; exit 1; }
  [ -f "$CELLS" ] || { echo "[0] missing $CELLS" >&2; exit 1; }

  # 200 reads off one locus is plenty - the two barcode spaces are disjoint,
  # so the answer is either ~all or exactly none.
  hit=$(samtools view "$SRC" chrX:100000000-100100000 2>/dev/null | head -200 \
        | grep -o 'CB:Z:[ACGT]*-1' | sed 's/CB:Z://' | sort -u \
        | grep -c -F -f <(cut -f1 "$CELLS") || true)
  echo "[0] $ID: $hit of the sampled CB tags are in $(basename "$CELLS")"
  if [ "$hit" -eq 0 ]; then
    echo "[0] $ID: ZERO overlap - the cells file is in the wrong barcode space." >&2
    echo "[0]     CB in this BAM is the CANONICAL (== gex) barcode, so use the" >&2
    echo "[0]     sinto_${ID}_rna*.txt files, NOT sinto_${ID}_atac*.txt." >&2
    exit 1
  fi
done

##### ============ stage A: filter the ATAC BAMs ============ #####
for ID in $SAMPLES; do
  SRC=$WORK/$ID/outs/atac_possorted_bam.bam
  FLT=$WORK/allelome_pseudobulk/bam/${ID}_atac_dedup.bam
  mkdir -p "$(dirname "$FLT")"
  if [ -s "$FLT.bai" ]; then
    echo "[A] $ID: filtered BAM present, skipping"
  else
    echo "[A] $(date +%T) $ID: filtering (-F 0x400 -q 30)"
    samtools view -b -F 0x400 -q 30 -@ 8 -o "$FLT" "$SRC" && samtools index -@ 8 "$FLT"
    echo "[A] $(date +%T) $ID: $(samtools view -c -@ 8 "$FLT") reads retained"
  fi
done

##### ---- helper: split with sinto, then score each BAM ---- #####
# $1 stage tag  $2 cells file  $3 tree name  $4 sample id
run_stage () {
  local TAG=$1 CELLS=$2 TREE=$3 ID=$4
  local FLT=$WORK/allelome_pseudobulk/bam/${ID}_atac_dedup.bam
  local BAMDIR=$SCR/$TREE/bam/$ID
  local OUTDIR=$SCR/$TREE/out/$ID

  [ -f "$CELLS" ] || { echo "[$TAG] $ID: missing $CELLS - skipping" >&2; return 0; }
  mkdir -p "$BAMDIR" "$OUTDIR"

  # DISTINCT GROUPS, not lines: the per-celltype file is one line per nucleus
  # with a group column, so wc -l would demand 1422 BAMs where sinto makes 7.
  local n_want; n_want=$(cut -f2 "$CELLS" | sort -u | wc -l)
  local n_have; n_have=$(ls "$BAMDIR"/*.bam 2>/dev/null | wc -l)
  if [ "$n_have" -lt "$n_want" ]; then
    echo "[$TAG] $(date +%T) $ID: sinto split into $n_want groups ($n_have BAMs already present)"
    conda activate sinto
    sinto filterbarcodes --bam "$FLT" --cells "$CELLS" --outdir "$BAMDIR" --nproc "$nproc"
    conda activate RNAseq
  else
    echo "[$TAG] $ID: $n_have BAMs already split, skipping sinto"
  fi

  local ilist=$SCR/$TREE/index_${ID}.txt; : > "$ilist"
  for b in "$BAMDIR"/*.bam; do
    [ -s "$b.bai" ] || echo "samtools index '$b'"
  done >> "$ilist"
  [ -s "$ilist" ] && parallel -j "$nproc" :::: "$ilist"

  # Allelome.PRO2 names its output <DATE>_<bambase>_<annot>_<minread>, so resume
  # tests for any such directory already holding a locus_table.
  local tlist=$SCR/$TREE/tasks_${ID}.txt; : > "$tlist"
  local done_n=0 todo_n=0
  for b in "$BAMDIR"/*.bam; do
    local base; base=$(basename "$b" .bam)
    if compgen -G "$OUTDIR/*_${base}_*/locus_table.txt" > /dev/null; then
      done_n=$((done_n+1))
    else
      echo "bash $PIPELINE -i '$b' -a $ANNOTATION -s $SNPS -r 1 -t 1 -o '$OUTDIR'" >> "$tlist"
      todo_n=$((todo_n+1))
    fi
  done
  echo "[$TAG] $ID: $done_n already scored, $todo_n to do"
  if [ "$todo_n" -gt 0 ]; then
    echo "[$TAG] $(date +%T) $ID: scoring"
    parallel -j "$nproc" :::: "$tlist"
    echo "[$TAG] $(date +%T) $ID: done, $(ls -d "$OUTDIR"/*/ 2>/dev/null | wc -l) output dirs"
  fi
}

##### ============ stage B: per cell type (fast, do first) ============ #####
# _rna_bycelltype, not _atac_bycelltype. See the barcode note in the header.
for ID in $SAMPLES; do
  run_stage B "$FIGS/sinto_${ID}_rna_bycelltype.txt" celltype "$ID"
done

##### ================= stage C: per nucleus ================= #####
for ID in $SAMPLES; do
  run_stage C "$QC/sinto_${ID}_rna.txt" pernucleus "$ID"
done

echo
echo "$(date +%T) all stages attempted. Consolidation and plotting are a"
echo "separate job's responsibility, so partial results still get plotted."
for TREE in celltype pernucleus; do
  for ID in $SAMPLES; do
    printf "  %-11s %-4s %6d locus tables\n" "$TREE" "$ID" \
      "$(find "$SCR/$TREE/out/$ID" -name locus_table.txt 2>/dev/null | wc -l)"
  done
done
