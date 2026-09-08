#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_allelome_pb
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
#SBATCH --time=8:00:00
#SBATCH --array=1-2

# Whole-sample pseudobulk Allelome.PRO2 on the multiome GEX BAMs. This is the
# FIRST allelic step and deliberately needs no cell typing: it answers whether
# the global chrX escape level reproduces the ~12.7% established by the
# snRNA-seq and spatial work, before any effort goes into WNN clustering and
# per-celltype sinto splits. If this does not reproduce, everything downstream
# is suspect.
#
# Array 1 = 9w, 2 = 78w.
#
# WHY chr_annotation_mm39.bed: it holds one row per chromosome (21 rows), so a
# single run yields both halves of the question at once -
#   chrX  -> the escape level, since CAST is the inactive X in every nucleus
#            (Xist is deleted on B6, so the B6 X cannot be silenced)
#   chr1..chr19 -> the mapping-bias baseline. Reads are aligned to a B6
#            reference, so CAST reads align slightly less well and every
#            allelic ratio is biased toward B6. Autosomes should sit at 0.5;
#            whatever they are off by is the correction chrX needs. Getting
#            that from the same run, same filters, same SNP set removes any
#            question of whether the control is comparable.
#
# WHY the _no_Xist SNP build: the B6 Xist allele is DELETED, so that interval
# carries no B6 sequence and scores as 100% CAST - a deletion artefact
# indistinguishable from perfect escape. Masking it is not optional here.
#
# WHY -F 0x400 -q 255, and what each actually does. These are NOT symmetric:
#
#   -q 255  changes the RESULT. STAR-unique only. mpileup applies no MAPQ
#           filter by default, so without this a multimapper overlapping a SNP
#           reaches the scoring and is assigned an allele outright.
#   -F 0x400  changes only the SPEED. Allelome.PRO2.sh line 238 runs
#           `samtools mpileup` with no --ff, so mpileup's default exclusion
#           mask UNMAP,SECONDARY,QCFAIL,DUP applies and duplicate-flagged reads
#           NEVER reach the allelic scoring. Pre-filtering them is therefore a
#           no-op for correctness. It is kept because these libraries are
#           98.58%/97.79% duplicate, so dropping them shrinks the BAM ~70x and
#           makes both the sinto split and the pileup dramatically faster.
#
# So Allelome.PRO2 has always been duplicate-excluded here, whatever an
# upstream filter says. total_reads in locus_table.txt is already a molecule-
# level count. That was measured the hard way - see the memory note
# allelome-pro2-drops-duplicates, where a 24h x 2 comparison produced
# byte-identical trees.
#
# NOTE for the ATAC counterpart, when it comes: -q 255 must NOT be reused.
# cellranger-arc aligns ATAC with BWA, whose MAPQ caps at 60, so -q 255 there
# returns zero reads WITHOUT erroring. Use -q 30.

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
GRCM39=$BASE/GRCm39
PIPELINE=$BASE/Allelome.LINK/00_src/Allelome.PRO2.sh
ANNOTATION=$GRCM39/chr_annotation_mm39.bed
SNPS=$GRCM39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed

case "${SLURM_ARRAY_TASK_ID:-0}" in
  1) ID=9w  ;;
  2) ID=78w ;;
  *) echo "unexpected SLURM_ARRAY_TASK_ID '${SLURM_ARRAY_TASK_ID:-}'" >&2; exit 1 ;;
esac

WORK=$BASE/adult_aged_multiome
SRC=$WORK/$ID/outs/gex_possorted_bam.bam
FLT_DIR=$WORK/allelome_pseudobulk/bam
OUT=$WORK/allelome_pseudobulk/$ID
FLT=$FLT_DIR/${ID}_gex_dedup.bam

for f in "$PIPELINE" "$ANNOTATION" "$SNPS" "$SRC"; do
  [ -f "$f" ] || { echo "Missing file: $f" >&2; exit 1; }
done
mkdir -p "$FLT_DIR" "$OUT"

##### ---- filter duplicates and multimappers ---- #####
# Resume on the INDEX, not the BAM: a job killed mid-`samtools view` leaves a
# truncated but non-empty BAM, and the && means it never got an index. The
# index is the only marker that says the pair of steps finished.
if [ -s "$FLT.bai" ]; then
  echo "[$(date +%T)] $ID: filtered BAM already present, skipping filter"
else
  echo "[$(date +%T)] $ID: filtering $SRC"
  samtools view -b -F 0x400 -q 255 -@ 8 -o "$FLT" "$SRC" && samtools index -@ 8 "$FLT"
fi

##### ---- report what the filter actually cost ---- #####
# Logged rather than assumed: the 98.58%/97.79% duplicate rates come from
# cellranger's own summary, and this is the independent check on them. It also
# gives the informative depth that per-nucleus power rests on.
echo "[$(date +%T)] $ID: read counts"
n_raw=$(samtools view -c -@ 8 "$SRC")
n_flt=$(samtools view -c -@ 8 "$FLT")
n_flt_x=$(samtools view -c -@ 8 "$FLT" chrX)
printf "   raw            %12d\n" "$n_raw"
printf "   after filter   %12d  (%.2f%% retained)\n" "$n_flt" \
  "$(awk -v a="$n_flt" -v b="$n_raw" 'BEGIN{print (b?100*a/b:0)}')"
printf "   after filter, chrX %8d\n" "$n_flt_x"

##### ---- Allelome.PRO2 ---- #####
# No -r/-t, matching Allelome.LINK.pseudobulk.slurm: pseudobulk uses the
# defaults rather than the -r 1 -t 1 that the single-cell runs need.
echo "[$(date +%T)] $ID: Allelome.PRO2"
bash "$PIPELINE" \
  -i "$FLT" \
  -a "$ANNOTATION" \
  -s "$SNPS" \
  -o "$OUT"

echo "[$(date +%T)] $ID: done. locus_table.txt under $OUT"
