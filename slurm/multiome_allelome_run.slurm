#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_allelome
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

# Allelome.PRO2 over the QC-passing multiome nuclei, both groupings, in one
# overnight pass. Stages run cheapest-first so there is a usable result even if
# the wall clock runs out:
#
#   A  filter both GEX BAMs once            (-F 0x400 -q 255)
#   B  PER CELL TYPE  - 14 BAMs, minutes    <- guaranteed to land
#   C  PER NUCLEUS    - 4547 BAMs, hours    <- the actual question
#
# EVERY STAGE IS RESUMABLE. Each unit is skipped when its output already
# exists, so a timeout plus a resubmit continues rather than restarting. That
# matters more than usual here: the per-nucleus cost was never measured
# (the pilot was written but never run), so stage C's duration is genuinely
# unknown and may need a second night.
#
# WHY chr_annotation_mm39.bed: one row per chromosome, so a single run yields
# chrX (the escape readout, since CAST is the inactive X in every nucleus) AND
# chr1-19 (the mapping-bias baseline - reads are aligned to a B6 reference, so
# CAST reads align slightly worse and every ratio is pulled toward B6).
# Autosomes should sit at 0.5; the offset is the correction chrX needs.
#
# WHY _no_Xist: the B6 Xist allele is DELETED, so that interval carries no B6
# sequence and would score as 100% CAST - a deletion artefact indistinguishable
# from perfect escape.
#
# WHY -F 0x400 -q 255: these libraries are 98.58%/97.79% duplicate, so read-
# level scoring would let one amplified molecule vote ~70 times. -q 255 is
# STAR-unique; a multimapper over a SNP is assigned an allele outright. NOTE
# this is the GEX BAM only - the ATAC BAM is BWA-aligned and would need -q 30.
#
# Per-nucleus BAMs and Allelome output live on SCRATCH at a STABLE path (no job
# id), so the plotting job can consolidate whatever finished even if this job
# was killed. Only small text reaches DSS.

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
GRCM39=$BASE/GRCm39
WORK=$BASE/adult_aged_multiome
PIPELINE=$BASE/Allelome.LINK/00_src/Allelome.PRO2.sh
ANNOTATION=$GRCM39/chr_annotation_mm39.bed
SNPS=$GRCM39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed
FIGS=$WORK/figures_signac        # where 03 wrote the per-celltype sinto files
QC=$WORK/qc                      # where 01 wrote the per-nucleus sinto files
SAMPLES="9w 78w"

for f in "$PIPELINE" "$ANNOTATION" "$SNPS"; do
  [ -f "$f" ] || { echo "Missing: $f" >&2; exit 1; }
done
if [ -z "${SCRATCH:-}" ]; then echo "SCRATCH is not set - aborting" >&2; exit 1; fi

SCR=$SCRATCH/multiome_allelome
mkdir -p "$SCR"
nproc=$((SLURM_CPUS_PER_TASK - 2)); [ "$nproc" -lt 1 ] && nproc=1
module load parallel
shopt -s nullglob
echo "scratch: $SCR   parallel jobs: $nproc"

##### ================= stage A: filter the GEX BAMs ================= #####
conda activate RNAseq
for ID in $SAMPLES; do
  SRC=$WORK/$ID/outs/gex_possorted_bam.bam
  FLT=$WORK/allelome_pseudobulk/bam/${ID}_gex_dedup.bam
  mkdir -p "$(dirname "$FLT")"
  if [ -s "$FLT.bai" ]; then
    echo "[A] $ID: filtered BAM present, skipping"
  else
    echo "[A] $(date +%T) $ID: filtering"
    samtools view -b -F 0x400 -q 255 -@ 8 -o "$FLT" "$SRC" && samtools index -@ 8 "$FLT"
    echo "[A] $(date +%T) $ID: $(samtools view -c -@ 8 "$FLT") reads retained"
  fi
done

##### ---- helper: split with sinto, then score each BAM ---- #####
# $1 stage tag  $2 cells file pattern dir  $3 cells filename  $4 tree name
run_stage () {
  local TAG=$1 CELLS=$2 TREE=$3 ID=$4
  local FLT=$WORK/allelome_pseudobulk/bam/${ID}_gex_dedup.bam
  local BAMDIR=$SCR/$TREE/bam/$ID
  local OUTDIR=$SCR/$TREE/out/$ID

  [ -f "$CELLS" ] || { echo "[$TAG] $ID: missing $CELLS - skipping" >&2; return 0; }
  mkdir -p "$BAMDIR" "$OUTDIR"

  local n_want; n_want=$(wc -l < "$CELLS")
  local n_have; n_have=$(ls "$BAMDIR"/*.bam 2>/dev/null | wc -l)
  if [ "$n_have" -lt "$n_want" ]; then
    echo "[$TAG] $(date +%T) $ID: sinto split, $n_want groups ($n_have already present)"
    conda activate sinto
    sinto filterbarcodes --bam "$FLT" --cells "$CELLS" --outdir "$BAMDIR" --nproc "$nproc"
    conda activate RNAseq
  else
    echo "[$TAG] $ID: $n_have BAMs already split, skipping sinto"
  fi

  # index anything missing an index
  local ilist=$SCR/$TREE/index_${ID}.txt; : > "$ilist"
  for b in "$BAMDIR"/*.bam; do
    [ -s "$b.bai" ] || echo "samtools index '$b'"
  done >> "$ilist"
  [ -s "$ilist" ] && parallel -j "$nproc" :::: "$ilist"

  # score anything without a locus_table. Allelome.PRO2 names its output
  # <DATE>_<bambase>_<annot>_<minread>, so resume tests for any such directory
  # that already holds a locus_table rather than for a fixed name.
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
for ID in $SAMPLES; do
  run_stage B "$FIGS/sinto_${ID}_rna_bycelltype.txt" celltype "$ID"
done

##### ================= stage C: per nucleus ================= #####
for ID in $SAMPLES; do
  run_stage C "$QC/sinto_${ID}_rna.txt" pernucleus "$ID"
done

echo
echo "$(date +%T) all stages attempted. Consolidation and plotting are the"
echo "plotting job's responsibility, so partial results still get plotted."
for TREE in celltype pernucleus; do
  for ID in $SAMPLES; do
    printf "  %-11s %-4s %6d locus tables\n" "$TREE" "$ID" \
      "$(find "$SCR/$TREE/out/$ID" -name locus_table.txt 2>/dev/null | wc -l)"
  done
done
