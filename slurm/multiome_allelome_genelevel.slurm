#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_allelome_gene
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=1-2
#
# Re-score the multiome PSEUDOBULK BAMs at GENE level instead of chromosome
# level. One array task per sample; GNU parallel over the whole-sample BAM plus
# the per-cell-type BAMs.
#
#   sbatch slurm/multiome_allelome_genelevel.slurm [allelome_tree]
#
# allelome_tree defaults to $SCRATCH/multiome_allelome and is where stage B of
# multiome_allelome_run.slurm left celltype/bam/<id>/*.bam. Pass the DSS copy
# once that exists. POSITIONAL, not an environment variable: every slurm header
# in this repo carries --export=NONE, so an exported variable never arrives.
#
# WHY THIS EXISTS. Every allelic run so far used chr_annotation_mm39.bed, one
# row per chromosome, so the entire allelic result was a single escape number
# per sample per cell type. That number is dominated by whichever chrX genes
# are most highly expressed and cannot distinguish "12.7% of chrX reads escape"
# from "three genes escape strongly and the rest are silent" - the same pooling
# artefact as pooling reads across cell types, one level up. Nothing downstream
# in the plan (escape vs non-escape genes, RNA/ATAC concordance, peak-to-gene
# linkage) is answerable without per-gene ratios.
#
# WHY annotation_us_mm39_gene_level.bed AND NOT ..._chrX.bed. 39881 gene bodies
# over 21 chromosomes rather than 3227 chrX intervals. Two reasons, both from
# pseudobulk_celltype_allelome.slurm, which this mirrors: a chrX-only bed cannot
# answer a question about chromosome ends, and the autosomes are the only honest
# null for chrX - same B6-reference mapping bias, same repeat content, no
# monoallelic biology. It also keeps the mapping-bias correction and the escape
# measurement in ONE table under identical filters.
#
# WHY -r 1 -t 1, where the whole-sample chr-level job passed no -r/-t at all.
# A gene enters locus_table.txt on a single read, so the depth cutoff lives in
# ONE place, the R analysis, where the same table can be re-thresholded. Read
# depth is the dominant confounder of an allelic ratio - a 2-read gene is 0 or 1
# by construction, never 0.5 - so it must stay tunable rather than baked into
# the output directories. At chromosome level that did not matter; at gene level
# it decides which genes exist.
#
# WHY ITS OWN OUTPUT TREE, allelome_gene_level/, rather than reusing
# allelome_pseudobulk/ or the scratch tree. multiome_allelome_run.slurm tests
# for prior output with the loose glob "*_${base}_*/locus_table.txt", which
# matches ANY annotation's directory. Re-scoring in place with a new annotation
# would make every unit look already-done and the job would silently do nothing.
# The resume check below is tight instead: date is the only wildcard.
#
# The layout is deliberately the one OCM/Allelome.PRO2_pseudobulk_celltype uses,
# i.e. <root>/<group>/<date>_<name>_<annotation>_<minread>/, so the existing
# gene-level analysis in OCM_heart/ap2_pseudobulk_common.R reads it unchanged
# via its AP2_ROOT and GROUPS overrides.
#
# Allelome.PRO2.sh rewrites the chr prefix of the SNP and annotation beds IN
# PLACE, which is a data race across parallel workers on shared reference files.
# The pre-flight asserts all three already agree on "chr", making those steps
# read-only no-ops. Do not remove it.

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq
module load parallel
shopt -s nullglob

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
GRCM39=$BASE/GRCm39
WORK=$BASE/adult_aged_multiome
PIPELINE=$BASE/Allelome.LINK/00_src/Allelome.PRO2.sh
ANNOTATION=$GRCM39/annotation_us_mm39_gene_level.bed
SNPS=$GRCM39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed
ANNOT_NAME=$(basename "$ANNOTATION")
MIN_READ=1
TOTAL_READ=1

TREE=${1:-$SCRATCH/multiome_allelome}

case "${SLURM_ARRAY_TASK_ID:-0}" in
  1) ID=9w  ;;
  2) ID=78w ;;
  *) echo "unexpected SLURM_ARRAY_TASK_ID '${SLURM_ARRAY_TASK_ID:-}'" >&2; exit 1 ;;
esac

CT_BAMDIR=$TREE/celltype/bam/$ID
WHOLE_BAM=$WORK/allelome_pseudobulk/bam/${ID}_gex_dedup.bam
FINAL=$WORK/allelome_gene_level/$ID
SNPCOUNT=$WORK/allelome_gene_level/gene_level_snp_count.bed

[ -n "$SCRATCH" ] || { echo "SCRATCH is not set - aborting" >&2; exit 1; }
for f in "$PIPELINE" "$ANNOTATION" "$SNPS"; do
  [ -s "$f" ] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done
mkdir -p "$FINAL"

##### ---- collect the pseudobulk BAMs ---- #####
# The whole-sample BAM goes FIRST in the task list. It is far larger than any
# single cell type, so parallel must start it immediately or it becomes the
# straggler that sets the job's wall time.
bams=()
[ -s "$WHOLE_BAM" ] && bams+=("$WHOLE_BAM")
if [ -d "$CT_BAMDIR" ]; then
  for b in "$CT_BAMDIR"/*.bam; do bams+=("$b"); done
else
  echo "WARNING: no per-cell-type BAMs at $CT_BAMDIR" >&2
  echo "         stage B of multiome_allelome_run.slurm writes them; pass the" >&2
  echo "         tree location as \$1 if they have been moved off scratch." >&2
fi
n_bams=${#bams[@]}
[ "$n_bams" -gt 0 ] || { echo "ERROR: no pseudobulk BAMs found for $ID" >&2; exit 1; }
echo "$ID: $n_bams pseudobulk BAMs (1 whole-sample + per cell type)"

##### ---- pre-flight: chr prefix must already agree (see header) ---- #####
bam_chr=$(samtools view -H "${bams[0]}" | awk '$1=="@SQ"{print $2; exit}' | sed 's/^SN://')
case "$bam_chr" in
  chr*) ;;
  *) echo "ERROR: BAM contig '$bam_chr' has no 'chr' prefix; Allelome.PRO2.sh would rewrite the shared SNP and annotation beds in place, racing across parallel workers. Aborting." >&2; exit 1 ;;
esac
head -n 1 "$SNPS"       | grep -q '^chr' || { echo "ERROR: $SNPS lacks the 'chr' prefix (see above)" >&2; exit 1; }
head -n 1 "$ANNOTATION" | grep -q '^chr' || { echo "ERROR: $ANNOTATION lacks the 'chr' prefix (see above)" >&2; exit 1; }
echo "chr-prefix pre-flight OK (BAM contig $bam_chr)"

##### ---- index anything missing one ---- #####
for b in "${bams[@]}"; do
  [ -s "$b.bai" ] || { echo "indexing $(basename "$b")"; samtools index "$b"; }
done

##### ---- one-time: informative SNPs per gene ---- #####
# Gene-level ratios need this and chromosome-level ones did not. SNP density
# varies along the chromosome, so a gene with few usable SNPs is measured worse
# rather than being more biallelic - a confounder for any positional trend, and
# chrX distal windows are exactly where that bites. Counted once, joined in R.
# Both array tasks reach this together: write to a per-task temp and mv, which
# is atomic on one filesystem, so the loser overwrites with an identical file
# instead of a reader seeing a half-written one. No -sorted: that needs the SNP
# bed and the annotation to agree on chromosome ORDER, which "sorted" in the
# filename does not promise, and bedtools errors out rather than falling back.
if [ ! -s "$SNPCOUNT" ]; then
  echo "counting informative SNPs per gene -> $SNPCOUNT"
  tmp_sc=$SCRATCH/gene_level_snp_count.${SLURM_ARRAY_TASK_ID:-0}.$$.bed
  tmp_an=$SCRATCH/annot_sorted.${SLURM_ARRAY_TASK_ID:-0}.$$.bed
  sort -k1,1 -k2,2n "$ANNOTATION" > "$tmp_an"
  if bedtools intersect -c -a "$tmp_an" -b "$SNPS" > "$tmp_sc"; then
    mv -f "$tmp_sc" "$SNPCOUNT"
  else
    echo "ERROR: bedtools intersect failed" >&2
    rm -f "$tmp_sc" "$tmp_an"; exit 1
  fi
  rm -f "$tmp_an"
fi

##### ---- run ---- #####
SCRDIR=$SCRATCH/multiome_gene_level_${ID}_${SLURM_JOBID}
rm -rf "$SCRDIR"; mkdir -p "$SCRDIR"
export TMPDIR="$SCRDIR"          # keep GNU sort's spill out of a shared /tmp

nproc=${SLURM_CPUS_PER_TASK:-8}
[ "$n_bams" -lt "$nproc" ] && nproc=$n_bams

# Resume on the TIGHT glob: date is the only wildcard. A loose *_${base}_*
# would match this sample's chromosome-level output and skip every unit.
TASKS=$FINAL/task-list.txt; : > "$TASKS"
n_done=0
for b in "${bams[@]}"; do
  base=$(basename "$b" .bam)
  if compgen -G "$FINAL"/*_"${base}"_"${ANNOT_NAME}"_"${MIN_READ}"/locus_table.txt > /dev/null; then
    n_done=$((n_done + 1))
  else
    echo "bash $PIPELINE -i '$b' -a $ANNOTATION -s $SNPS -r $MIN_READ -t $TOTAL_READ -o '$SCRDIR'" >> "$TASKS"
  fi
done
n_todo=$(wc -l < "$TASKS")
echo "$ID: $n_done already scored, $n_todo to do, parallel -j $nproc"

if [ "$n_todo" -gt 0 ]; then
  # No --halt: a very small cell-type pseudobulk can legitimately fail with an
  # empty pileup, and that must not take the others down with it.
  parallel --joblog "$FINAL/joblog.txt" -j "$nproc" :::: "$TASKS"
fi

##### ---- harvest the small outputs, drop the trimmed BAMs and pileups ---- #####
for d in "$SCRDIR"/*/; do
  [ -d "$d" ] || continue
  dest=$FINAL/$(basename "$d"); mkdir -p "$dest"
  for f in locus_table.txt read_count_per_SNP.txt; do
    [ -s "$d$f" ] && mv "$d$f" "$dest/"
  done
  for f in "$d"*.log;           do [ -s "$f" ] && mv "$f" "$dest/"; done
  for f in "$d"BED_files/*.bed; do [ -s "$f" ] && mv "$f" "$dest/"; done
done

##### ---- judge success by OUTPUT, not exit status ---- #####
# Allelome.PRO2.sh answers its own argument checks with `exit 0` and has no
# set -e, so its exit code cannot report a failure. A run counts as successful
# only if it left a locus_table.txt behind.
n_fail=0; failed=""
for b in "${bams[@]}"; do
  base=$(basename "$b" .bam)
  compgen -G "$FINAL"/*_"${base}"_"${ANNOT_NAME}"_"${MIN_READ}"/locus_table.txt > /dev/null \
    || { n_fail=$((n_fail + 1)); failed="$failed $base"; }
done
if [ "$n_fail" -gt 0 ]; then
  echo "WARNING: $n_fail of $n_bams pseudobulks produced no locus_table.txt:" >&2
  for f in $failed; do echo "  $f" >&2; done
  echo "  A cell type with very few nuclei has almost no SNP-overlapping reads" >&2
  echo "  and cannot be scored; check its .log under $SCRDIR." >&2
fi
if [ "$n_fail" -eq "$n_bams" ]; then
  echo "ERROR: every Allelome.PRO2 run failed for $ID" >&2
  exit 1
fi

##### ---- drop the scratch intermediates ---- #####
# Allelome.PRO2 leaves a trimmed BAM and a pileup behind for EVERY unit it
# scores, and they are close to a fixed size per run regardless of how small the
# input BAM is. The per-nucleus stage of multiome_allelome_run.slurm has no such
# cleanup, and its output tree reached 6.7 TB to yield locus tables totalling a
# few tens of MB. Sixteen pseudobulk units would be far smaller than that, but
# the harvest above has already taken everything worth keeping, so there is no
# reason to leave them. Guarded on the harvest having produced something.
if [ "$n_fail" -lt "$n_bams" ]; then
  rm -rf "$SCRDIR"
fi

echo "$ID done -> $FINAL"
