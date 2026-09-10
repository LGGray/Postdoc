#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J pseudobulk_celltype_allelome
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16   # serial_std range is 1-16; sbatch REJECTS 17+
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --array=1-4
#
# Allelome.PRO2 on the per-cell-type pseudobulk BAMs, one array task per
# cohort (9w / 78w / Sham / TAC), GNU parallel over the ~11 cell types.
#
# serial_std, not cm4_tiny. A cohort has at most 11 cell types and every
# Allelome.PRO2 run is single-threaded, so the job needs -j 11 and no more:
# serial_std's 1-16 core range is a better fit than cm4_tiny's 17-112, where
# 32 was requested and half of it sat idle. serial_std also places
# essentially immediately, whereas the cm4_tiny attempt sat PD (Priority) --
# the same trade documented in submit_sinto_tiles_chain.sh. 150G fits well
# inside serial_std's 1000 GiB node limit, and 4 array tasks x 16 cores = 64
# is under the 96-core sum-over-running-jobs cap.
#
#   sbatch slurm/pseudobulk_celltype_allelome.slurm
#
# Prerequisite: slurm/pseudobulk_celltype_bams.slurm has written
# pseudobulk_celltype/<cohort>/*.bam + bamfiles.txt.
#
# Two deliberate differences from the per-cell runs in scAllelome_*.slurm:
#
#  * annotation_us_mm39_gene_level.bed, not annotation_us_mm39_chrX.bed.
#    39881 gene bodies over 21 chromosomes. The chrX bed cannot answer a
#    question about chromosome ends, and the autosomes are also the only
#    honest null for it: they carry the same B6-reference mapping bias and
#    the same subtelomeric repeat content as chrX, with no monoallelic
#    biology.
#
#  * -t 1, so a gene enters locus_table.txt on a single read. Pseudobulk has
#    real depth, so this is not the necessity it was per cell -- it is so the
#    depth cutoff lives in ONE place, the R analysis, where the same table can
#    be re-thresholded. Read depth is the dominant confounder of an allelic
#    ratio (a 2-read gene is 0 or 1 by construction, never 0.5), so it must
#    stay a tunable, not something baked into 44 output directories.
#
# The SNP bed is the no-Xist build, per the house rule: these animals carry the
# Xist deletion on the B6 X, so SNPs inside the Xist span cannot be read as
# allelic imbalance.
#
# Allelome.PRO2.sh rewrites the chr prefix of the SNP and annotation files IN
# PLACE, which would be a data race across 11 parallel workers on shared
# reference files. The pre-flight below asserts all three already agree on the
# "chr" prefix, which makes those steps read-only no-ops. Do not remove it.

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq
module load parallel

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM

case "${SLURM_ARRAY_TASK_ID:-0}" in
	1) cohort="9w"   ;;
	2) cohort="78w"  ;;
	3) cohort="Sham" ;;
	4) cohort="TAC"  ;;
	*) echo "ERROR: unknown SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2; exit 1 ;;
esac

root="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2"
pipeline="$root/Allelome.LINK/00_src/Allelome.PRO2.sh"
annotation="$root/GRCm39/annotation_us_mm39_gene_level.bed"
snps="$root/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed"

annotation_name=$(basename "$annotation")
bam_dir="pseudobulk_celltype/${cohort}"
input="$bam_dir/bamfiles.txt"
final_dir="$root/OCM/Allelome.PRO2_pseudobulk_celltype/${cohort}"
min_read=1
total_read=1

if [ -z "$SCRATCH" ]; then echo "SCRATCH is not set - aborting" >&2; exit 1; fi
for f in "$pipeline" "$annotation" "$snps" "$input"; do
	[ -s "$f" ] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done

# ---- pre-flight: chr prefix must already agree (see header) ----------------
one_bam="$bam_dir/$(head -n 1 "$input")"
[ -s "$one_bam" ] || { echo "ERROR: missing $one_bam" >&2; exit 1; }
bam_chr=$(samtools view -H "$one_bam" | awk '$1=="@SQ"{print $2; exit}' | sed 's/^SN://')
case "$bam_chr" in chr*) ;; *) echo "ERROR: BAM contig '$bam_chr' has no 'chr' prefix; Allelome.PRO2.sh would rewrite the shared SNP and annotation beds in place, racing across parallel workers. Aborting." >&2; exit 1 ;; esac
head -n 1 "$snps"       | grep -q '^chr' || { echo "ERROR: $snps lacks the 'chr' prefix (see above)" >&2; exit 1; }
head -n 1 "$annotation" | grep -q '^chr' || { echo "ERROR: $annotation lacks the 'chr' prefix (see above)" >&2; exit 1; }
echo "chr-prefix pre-flight OK (BAM contig $bam_chr)"

# ---- one-time: SNPs per gene, the other confounder of a positional trend ---
# Informative SNP density varies along a chromosome, and subtelomeric and
# pericentromeric sequence is repeat-rich. A gene with few usable SNPs is
# measured worse, not necessarily biallelically. Counted once, joined in R.
# Written to a per-task temp then mv'd into place: all four array tasks reach
# this at once, and mv on the same filesystem is atomic, so the losers
# overwrite with an identical file instead of a reader seeing a half-written
# one. No -sorted: that would require the 630 MB SNP bed and the annotation to
# agree on chromosome ORDER, which "sorted" in its filename does not promise
# (chr1, chr2, ... chr19, chrX vs sort's chr1, chr10, chr11, ...), and
# bedtools errors out rather than falling back. Loading the SNPs unsorted
# costs a few GB of the 150 G once.
snpcount="pseudobulk_celltype/gene_level_snp_count.bed"
if [ ! -s "$snpcount" ]; then
	echo "counting informative SNPs per gene -> $snpcount"
	tmp_snpcount="$SCRATCH/gene_level_snp_count.${SLURM_ARRAY_TASK_ID:-0}.$$.bed"
	sort -k1,1 -k2,2n "$annotation" > "$SCRATCH/annot_sorted_$$.bed"
	if bedtools intersect -c -a "$SCRATCH/annot_sorted_$$.bed" -b "$snps" > "$tmp_snpcount"; then
		mv -f "$tmp_snpcount" "$snpcount"
	else
		echo "ERROR: bedtools intersect failed" >&2
		rm -f "$tmp_snpcount" "$SCRATCH/annot_sorted_$$.bed"
		exit 1
	fi
	rm -f "$SCRATCH/annot_sorted_$$.bed"
fi

# ---- run ------------------------------------------------------------------
scratch_dir="$SCRATCH/pseudobulk_allelome_${cohort}_${SLURM_JOBID}"
rm -rf "$scratch_dir"
mkdir -p "$scratch_dir"
mkdir -p "$final_dir"

# keep GNU sort's spill inside the job's scratch, not a shared /tmp
export TMPDIR="$scratch_dir"

n_bams=$(wc -l < "$input")
nproc="${SLURM_CPUS_PER_TASK:-8}"
if [ "$n_bams" -lt "$nproc" ]; then nproc="$n_bams"; fi
echo "cohort=$cohort  bams=$n_bams  parallel -j $nproc"

tasks="$final_dir/task-list.txt"
while IFS= read -r line; do
	echo "bash $pipeline -i $bam_dir/$line -a $annotation -s $snps -r $min_read -t $total_read -o $scratch_dir"
done < "$input" > "$tasks"

# No --halt: the very small pseudobulks (78w stressed CM is n=1 nucleus) can
# legitimately fail with an empty pileup, and that must not take the other ten
# cell types down with it. Failures are collected from the joblog and reported;
# the job only fails outright if EVERY pseudobulk failed.
joblog="$final_dir/joblog.txt"
parallel --joblog "$joblog" -j "$nproc" :::: "$tasks"

# Success is judged by OUTPUT, not by exit status. Allelome.PRO2.sh answers
# its own argument checks with `exit 0` ("Error: Sample bam file ... can not be
# found." then exit 0) and has no set -e, so its exit code cannot be trusted to
# report a failure. A run counts as successful only if it left a
# locus_table.txt behind. The joblog is still written, for wall times.
n_run=0; n_fail=0; failed=""
while IFS= read -r line; do
	n_run=$((n_run + 1))
	sample_n="${line%.bam}"
	# Allelome.PRO2.sh names its job dir <date>_<sample>_<annotation>_<minreads>,
	# so only the date is a wildcard here -- a loose *_${sample_n}_* could match
	# a different cell type's directory and hide a failure.
	if ! ls "$scratch_dir"/*_"${sample_n}"_"$annotation_name"_"$min_read"/locus_table.txt >/dev/null 2>&1; then
		n_fail=$((n_fail + 1))
		failed="$failed $sample_n"
	fi
done < "$input"

if [ "$n_fail" -gt 0 ]; then
	echo "WARNING: $n_fail of $n_run pseudobulks produced no locus_table.txt:" >&2
	for f in $failed; do echo "  $f" >&2; done
	echo "  Check that pseudobulk's .log under $scratch_dir; one with almost no" >&2
	echo "  nuclei has almost no SNP-overlapping reads and cannot be scored." >&2
fi
if [ "$n_run" -gt 0 ] && [ "$n_fail" -eq "$n_run" ]; then
	echo "ERROR: every Allelome.PRO2 run failed for $cohort" >&2
	rc=1
else
	rc=0
fi

# ---- keep the small outputs, drop the trimmed BAMs and pileups ------------
for sample_dir in "$scratch_dir"/*/; do
	[ -d "$sample_dir" ] || continue
	sample_name=$(basename "$sample_dir")
	dest_dir="$final_dir/$sample_name"
	mkdir -p "$dest_dir"
	for f in locus_table.txt read_count_per_SNP.txt; do
		[ -s "$sample_dir$f" ] && mv "$sample_dir$f" "$dest_dir/"
	done
	# the .log carries the read/SNP accounting AND is named <celltype>.log,
	# which is how the R analysis recovers the cell type from the job dir
	for f in "$sample_dir"*.log; do [ -s "$f" ] && mv "$f" "$dest_dir/"; done
	for f in "$sample_dir"BED_files/*.bed; do [ -s "$f" ] && mv "$f" "$dest_dir/"; done
done

echo "$cohort done -> $final_dir"
exit "$rc"
