#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J pseudobulk_celltype_bams
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
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#
# Build one pseudobulk BAM per (cohort, cell type) from the OCM hearts, with
# scDblFinder doublets excluded, ready for Allelome.PRO2.
#
# --export=NONE, so the mode is POSITIONAL, not an environment variable:
#   sbatch slurm/pseudobulk_celltype_bams.slurm merge     <- default, recommended
#   sbatch slurm/pseudobulk_celltype_bams.slurm sinto
#   sbatch slurm/pseudobulk_celltype_bams.slurm merge 200   <- smaller merge batch
#
# serial_std, not cm4_tiny: samtools merge is I/O bound, so 16 threads and 18
# are indistinguishable here, and every stage of this chain has to live on the
# same cluster -- an LRZ --dependency does not cross clusters.
#
# Prerequisite (login node, seconds):
#   cd $OCM && Rscript $POSTDOC_ROOT/OCM_heart/pseudobulk_celltype_ids.R
#
# Two routes to the same reads:
#
# merge  Concatenates the per-cell BAMs that sinto filterbarcodes already
#        produced in sinto/<cohort>/, taking only the singlet barcodes. Those
#        per-cell sets are an exact 1:1 match to the QC'd cell list
#        (2376 / 2489 / 6927 / 8093 BAMs vs the same counts in
#        scDblFinder_per_cell.txt), and sinto assigns each read to exactly one
#        barcode, so merging the singlets of one cell type gives byte-for-byte
#        the same read set as re-running filterbarcodes with a cell type map.
#        This is the default because it needs no other input.
#
# sinto  Runs sinto filterbarcodes directly on the cellranger per-sample BAM
#        with the cell type map. Equivalent, but Hearts_OCM/ is currently only
#        present as Hearts_OCM.tar.gz (79 GB), so this route needs the archive
#        unpacked first. The script checks and tells you.
#
# samtools merge in BATCHES, not one merge over the whole group: merge holds
# every input file open at once and the largest group here is ~2000 per-cell
# BAMs, past the usual 1024 open-file limit. Batching keeps it well under, and
# because the per-cell BAMs are already coordinate-sorted, merge preserves that
# order -- so the result needs no sort pass and Allelome.PRO2 can be given
# -so 1 (its default). samtools cat would be one pass instead of two, but it
# requires byte-identical headers across every input and produces an unsorted
# BAM that would then have to be sorted anyway.

module load slurm_setup

MODE="${1:-merge}"
if [ "$MODE" != "merge" ] && [ "$MODE" != "sinto" ]; then
	echo "ERROR: mode must be 'merge' or 'sinto', got '$MODE'" >&2
	exit 1
fi

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
if [ "$MODE" = "sinto" ]; then
	conda activate sinto
else
	conda activate RNAseq
fi

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM

case "${SLURM_ARRAY_TASK_ID:-0}" in
	1) cohort="9w"   ; cr_bam="Hearts_OCM/outs/per_sample_outs/He_9w_snRNA_XBxC_het_XX/count/sample_alignments.bam" ;;
	2) cohort="78w"  ; cr_bam="Hearts_OCM/outs/per_sample_outs/He_78w_snRNA_XBxC_het_XX/count/sample_alignments.bam" ;;
	3) cohort="Sham" ; cr_bam="Hearts_OCM/outs/per_sample_outs/He_Sham_28d_snRNA_XBxC_het_XX/count/sample_alignments.bam" ;;
	4) cohort="TAC"  ; cr_bam="Hearts_OCM/outs/per_sample_outs/He_TAC_28d_snRNA_XBxC_het_XX/count/sample_alignments.bam" ;;
	*) echo "ERROR: unknown SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2; exit 1 ;;
esac

map="pseudobulk_celltype/sinto_celltype_${cohort}.txt"
out_dir="pseudobulk_celltype/${cohort}"
nproc="${SLURM_CPUS_PER_TASK:-8}"
# per-cell BAMs merged at once; stays clear of the open-file limit (ulimit -n)
# --export=NONE, so this is the second POSITIONAL argument, not an env var
BATCH="${2:-400}"

echo "cohort=$cohort  mode=$MODE  nproc=$nproc  merge batch=$BATCH"
ulimit -n 2>/dev/null | sed 's/^/open-file limit: /'  || true

if [ ! -s "$map" ]; then
	echo "ERROR: $map missing. Run OCM_heart/pseudobulk_celltype_ids.R first." >&2
	exit 1
fi
if [ -z "$SCRATCH" ]; then echo "SCRATCH is not set - aborting" >&2; exit 1; fi

scratch_dir="$SCRATCH/pseudobulk_celltype_${cohort}_${SLURM_JOBID}"
rm -rf "$scratch_dir"
mkdir -p "$scratch_dir" "$out_dir"

# groups to a file and read line by line rather than `for g in $groups`:
# word splitting is one shell setting away from breaking, and this loop is the
# thing that decides which nuclei end up in which pseudobulk BAM.
group_file="$scratch_dir/groups.txt"
cut -f2 "$map" | sort -u > "$group_file"
echo "groups in $cohort:"; sed 's/^/  /' "$group_file"

if [ "$MODE" = "sinto" ]; then
	if [ ! -s "$cr_bam" ]; then
		echo "ERROR: $cr_bam not found." >&2
		echo "  Hearts_OCM/ is archived as Hearts_OCM.tar.gz; unpack it, or use the" >&2
		echo "  merge route:  sbatch slurm/pseudobulk_celltype_bams.slurm merge" >&2
		exit 1
	fi
	export OMP_NUM_THREADS="$nproc"
	sinto filterbarcodes --bam "$cr_bam" --cells "$map" \
		--outdir "$scratch_dir" --nproc "$nproc" || exit 1

	# sinto writes coordinate-sorted output already; just move and index.
	while IFS= read -r g <&3; do
		if [ ! -s "$scratch_dir/${g}.bam" ]; then
			echo "ERROR: sinto produced no BAM for group $g" >&2; exit 1
		fi
		mv "$scratch_dir/${g}.bam" "$out_dir/${g}.bam" || exit 1
		samtools index -@ "$nproc" "$out_dir/${g}.bam" || exit 1
		echo "  $g: $(samtools view -c -@ "$nproc" "$out_dir/${g}.bam") reads"
	done 3< "$group_file"
else
	while IFS= read -r g <&3; do
		list="$scratch_dir/${g}.bamlist"
		# barcodes of this group -> the per-cell BAM sinto already wrote
		awk -F'\t' -v g="$g" -v co="$cohort" \
			'$2 == g { print "sinto/" co "/" co "_" $1 ".bam" }' "$map" > "$list"

		n_in=$(wc -l < "$list")
		if [ "$n_in" -eq 0 ]; then
			echo "ERROR: group $g has no barcodes in $map" >&2; exit 1
		fi

		# fail loudly on a missing per-cell BAM rather than silently
		# pseudobulking a subset of the cell type
		missing=$(while IFS= read -r f; do [ -s "$f" ] || echo "$f"; done < "$list")
		if [ -n "$missing" ]; then
			echo "ERROR: $g: per-cell BAMs missing from sinto/$cohort/:" >&2
			echo "$missing" | head -5 >&2
			echo "  ($(echo "$missing" | wc -l) of $n_in)" >&2
			exit 1
		fi

		echo "  $g: merging $n_in per-cell BAMs"
		if [ "$n_in" -le "$BATCH" ]; then
			samtools merge -f -@ "$nproc" -b "$list" "$out_dir/${g}.bam" || exit 1
		else
			batch_dir="$scratch_dir/${g}.batches"
			rm -rf "$batch_dir"; mkdir -p "$batch_dir"
			split -l "$BATCH" -d -a 4 "$list" "$batch_dir/part_"
			for part in "$batch_dir"/part_*; do
				samtools merge -f -@ "$nproc" -b "$part" "${part}.bam" || exit 1
			done
			ls "$batch_dir"/part_*.bam > "$batch_dir/level2.txt"
			echo "    $g: $(wc -l < "$batch_dir/level2.txt") batch BAMs -> final"
			samtools merge -f -@ "$nproc" -b "$batch_dir/level2.txt" \
				"$out_dir/${g}.bam" || exit 1
			rm -rf "$batch_dir"
		fi
		samtools index -@ "$nproc" "$out_dir/${g}.bam" || exit 1
		echo "  $g: $n_in nuclei, $(samtools view -c -@ "$nproc" "$out_dir/${g}.bam") reads"
	done 3< "$group_file"
fi

# bamfiles.txt: basenames, one per line -- the same convention
# sinto/<cohort>/bamfiles.txt uses, and what the Allelome.PRO2 driver reads.
( cd "$out_dir" && ls *.bam | sort > bamfiles.txt )
echo "wrote $out_dir/bamfiles.txt ($(wc -l < "$out_dir/bamfiles.txt") BAMs)"

rm -rf "$scratch_dir"
echo "$cohort done"
