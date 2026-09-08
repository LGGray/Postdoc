#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_pilot
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
#SBATCH --time=8:00:00

# Per-nucleus Allelome.PRO2 PILOT. Splits a small number of QC-passing nuclei
# out of the GEX BAM and scores them, purely to measure the per-nucleus cost
# before committing to all 4547 (1422 at 9w + 3125 at 78w).
#
#   sbatch slurm/multiome_sinto_pilot.slurm            # 50 nuclei of 9w
#   sbatch slurm/multiome_sinto_pilot.slurm 100 78w    # 100 nuclei of 78w
#
# Measuring first rather than guessing is the same discipline
# scAllelome_dedup.slurm used ("Measured on 50 9w cells (1.14M reads)"), and it
# matters more here: the OCM per-cell runs were sized for cm4_tiny at 110 cpus,
# and serial_std gives 16.
#
# ONE DELIBERATE DEPARTURE FROM THE OCM ROUTE. scAllelome_dedup.slurm split the
# BAM first and then filtered each per-cell BAM, i.e. one samtools pass per
# nucleus. Here the whole-BAM filter happens ONCE, before the split. Both
# filters are per-read and barcode-independent, so the result is identical, but
# this is one pass over 23G instead of ~1400 passes over small files.
#
# Note what that filter is and is not doing. -q 255 changes the result. -F
# 0x400 does not: Allelome.PRO2.sh line 238 runs mpileup with no --ff, so
# mpileup's default UNMAP,SECONDARY,QCFAIL,DUP mask already excludes duplicates
# before scoring. It is kept purely for speed - at ~98% duplicate it shrinks
# the BAM about 70x. The
# filtered BAM is shared with multiome_allelome_pseudobulk.slurm, which writes
# it to the same path, so whichever job runs first pays for it.
#
# Per-cell BAMs go to SCRATCH, not DSS. 4547 nuclei means ~9000 small files
# with indexes, which is exactly what a shared parallel filesystem is worst at.
# Only the Allelome.PRO2 output, which is small text, is kept.

module load slurm_setup

PILOT_N="${1:-50}"
ID="${2:-9w}"
case "$ID" in 9w|78w) ;; *) echo "ID must be 9w or 78w, got '$ID'" >&2; exit 1 ;; esac

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
GRCM39=$BASE/GRCm39
WORK=$BASE/adult_aged_multiome
PIPELINE=$BASE/Allelome.LINK/00_src/Allelome.PRO2.sh
ANNOTATION=$GRCM39/chr_annotation_mm39.bed
SNPS=$GRCM39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed

SRC=$WORK/$ID/outs/gex_possorted_bam.bam
CELLS=$WORK/qc/sinto_${ID}_rna.txt
FLT=$WORK/allelome_pseudobulk/bam/${ID}_gex_dedup.bam
OUT=$WORK/allelome_pernucleus_pilot/$ID

for f in "$PIPELINE" "$ANNOTATION" "$SNPS" "$SRC" "$CELLS"; do
  [ -f "$f" ] || { echo "Missing file: $f" >&2; exit 1; }
done
if [ -z "${SCRATCH:-}" ]; then echo "SCRATCH is not set - aborting" >&2; exit 1; fi

mkdir -p "$(dirname "$FLT")" "$OUT"
SPLIT=$SCRATCH/multiome_pilot_${SLURM_JOB_ID}
mkdir -p "$SPLIT"
echo "scratch: $SPLIT"

nproc=$((SLURM_CPUS_PER_TASK - 2)); [ "$nproc" -lt 1 ] && nproc=1
t0=$(date +%s)

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

##### ---- stage 1: filter the whole BAM once ---- #####
# Resume on the INDEX, not the BAM: a job killed mid-view leaves a truncated
# but non-empty BAM, and the && means it never got an index.
conda activate RNAseq
if [ -s "$FLT.bai" ]; then
  echo "[stage 1] filtered BAM already present, skipping"
else
  echo "[stage 1] samtools view -F 0x400 -q 255 over $SRC"
  samtools view -b -F 0x400 -q 255 -@ 8 -o "$FLT" "$SRC" && samtools index -@ 8 "$FLT"
fi
t1=$(date +%s)

##### ---- stage 2: split the pilot nuclei out ---- #####
head -n "$PILOT_N" "$CELLS" > "$SPLIT/cells.txt"
n_cells=$(wc -l < "$SPLIT/cells.txt")
echo "[stage 2] sinto filterbarcodes for $n_cells nuclei"
conda activate sinto
sinto filterbarcodes --bam "$FLT" --cells "$SPLIT/cells.txt" \
  --outdir "$SPLIT/bam" --nproc "$nproc"
t2=$(date +%s)

##### ---- stage 3: index, then score ---- #####
conda activate RNAseq
module load parallel
shopt -s nullglob
ls "$SPLIT"/bam/*.bam > "$SPLIT/bamlist.txt" || true
n_bam=$(wc -l < "$SPLIT/bamlist.txt")
echo "[stage 3] $n_bam per-nucleus BAMs produced"
if [ "$n_bam" -eq 0 ]; then
  echo "ERROR: sinto produced no BAMs. The usual cause is a barcode-space" >&2
  echo "       mismatch - the GEX BAM is in gex_barcode space, the ATAC BAM" >&2
  echo "       in atac_barcode space, and they never coincide." >&2
  exit 1
fi
awk '{print "samtools index \x27" $0 "\x27"}' "$SPLIT/bamlist.txt" > "$SPLIT/index-list.txt"
parallel -j "$nproc" :::: "$SPLIT/index-list.txt"
t3=$(date +%s)

# -r 1 -t 1 as the OCM single-cell runs use: a nucleus carries too few reads
# for the pseudobulk defaults to admit anything.
awk -v p="$PIPELINE" -v a="$ANNOTATION" -v s="$SNPS" -v o="$SPLIT/allelome" \
  '{print "bash " p " -i \x27" $0 "\x27 -a " a " -s " s " -r 1 -t 1 -o " o}' \
  "$SPLIT/bamlist.txt" > "$SPLIT/task-list.txt"
mkdir -p "$SPLIT/allelome"
echo "[stage 4] Allelome.PRO2 over $n_bam nuclei, -j $nproc"
parallel -j "$nproc" :::: "$SPLIT/task-list.txt"
t4=$(date +%s)

##### ---- keep the small text output, report the cost ---- #####
n_out=$(find "$SPLIT/allelome" -name locus_table.txt | wc -l)
cp "$SPLIT/cells.txt" "$SPLIT/task-list.txt" "$OUT/" 2>/dev/null
tar -czf "$OUT/locus_tables_${n_cells}nuclei.tar.gz" -C "$SPLIT" allelome 2>/dev/null

echo
echo "======================= pilot cost ======================="
printf "  nuclei requested          %d\n" "$n_cells"
printf "  per-nucleus BAMs          %d\n" "$n_bam"
printf "  locus_table.txt produced  %d\n" "$n_out"
printf "  stage 1 filter whole BAM  %5d s  (once per sample, already cached now)\n" "$((t1-t0))"
printf "  stage 2 sinto split       %5d s\n" "$((t2-t1))"
printf "  stage 3 index             %5d s\n" "$((t3-t2))"
printf "  stage 4 Allelome.PRO2     %5d s\n" "$((t4-t3))"
split_score=$((t4-t1))
echo "  ---"
printf "  split+score, %d nuclei     %5d s  => %.2f s/nucleus at -j %d\n" \
  "$n_cells" "$split_score" \
  "$(awk -v a="$split_score" -v b="$n_cells" 'BEGIN{print (b?a/b:0)}')" "$nproc"
for total in 1422 3125 4547; do
  printf "  projected for %4d nuclei  %5.1f h\n" "$total" \
    "$(awk -v s="$split_score" -v n="$n_cells" -v t="$total" 'BEGIN{print (n? s/n*t/3600 :0)}')"
done
echo "=========================================================="
echo "NOTE: stage 2 does not scale linearly - sinto makes one pass over the"
echo "      BAM regardless of how many barcodes it is splitting, so the full"
echo "      run's split will cost far less per nucleus than this pilot's."
echo "      Stage 4 is the one that scales, so weigh it most."
