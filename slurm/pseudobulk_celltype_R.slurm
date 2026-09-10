#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J pseudobulk_celltype_R
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4    # serial_std range is 1-16, no minimum
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#
# The two R stages of the per-cell-type pseudobulk chain, so the whole thing
# can be one dependency chain instead of two sbatch calls with two logins in
# between. See OCM_heart/PSEUDOBULK_CELLTYPE.md.
#
# --export=NONE, so everything is POSITIONAL:
#   sbatch slurm/pseudobulk_celltype_R.slurm ids
#   sbatch slurm/pseudobulk_celltype_R.slurm analysis
#   sbatch slurm/pseudobulk_celltype_R.slurm analysis Allelic_ratio_results_nodoublet 10
#
#   $1  stage          ids | analysis
#   $2  RESULTS_ROOT   default Allelic_ratio_results
#   $3  TERMINAL_MB    analysis only, default 5
#   $4  N_PERM         analysis only, default 2000
#
# serial_std. Both stages are single-threaded, and serial_std's 1-16 core
# range has NO minimum, so a small job costs 4 cores rather than the 17 that
# cm4_tiny's floor forces. Every stage of this chain must sit on the same
# cluster in any case: an LRZ --dependency does not cross clusters.
#
# If you ever move this back to cm4_tiny, --cpus-per-task must be >= 17.
# That partition's documented range is 17-112 PHYSICAL cores and sbatch
# rejects anything under it with "QOSMinCpuNotSatisfied" rather than rounding
# up. The partition table is the authority, not a neighbouring script:
# https://doku.lrz.de/job-processing-on-the-linux-cluster-10745970.html
#
# ids is seconds and needs only base R; analysis has to hold ~1.8 M gene rows
# plus the permutation draws, hence the 48G, and it sources
# allelic_ratio/00_functions.R, which is why this is seurat_env and not the
# repo-wide RNAseq env.

module load slurm_setup

stage="${1:-}"
results_root="${2:-Allelic_ratio_results}"
terminal_mb="${3:-5}"
n_perm="${4:-2000}"

case "$stage" in
	ids|analysis) ;;
	*) echo "ERROR: stage must be 'ids' or 'analysis', got '$stage'" >&2; exit 1 ;;
esac

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate seurat_env

POSTDOC_ROOT="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc"
export POSTDOC_ROOT
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM

echo "stage=$stage  RESULTS_ROOT=$results_root"

if [ "$stage" = "ids" ]; then
	Rscript "$POSTDOC_ROOT/OCM_heart/pseudobulk_celltype_ids.R" || exit 1

	# The chain's next stage reads these, so fail here rather than letting the
	# BAM build start and discover it group by group.
	for co in 9w 78w Sham TAC; do
		f="pseudobulk_celltype/sinto_celltype_${co}.txt"
		[ -s "$f" ] || { echo "ERROR: $f was not written" >&2; exit 1; }
		echo "  $f: $(wc -l < "$f") singlet barcodes, $(cut -f2 "$f" | sort -u | wc -l) groups"
	done
	exit 0
fi

# ---- analysis --------------------------------------------------------------
# Pre-flight the one package that is easy to be missing, so a 3-hour
# allocation is not spent discovering it 30 seconds in.
Rscript -e 'if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("tidyr is not installed in seurat_env; 13_pseudobulk_celltype_ar.R needs it")
}' || exit 1

export RESULTS_ROOT="$results_root"
export TERMINAL_MB="$terminal_mb"
export N_PERM="$n_perm"
echo "TERMINAL_MB=$terminal_mb  N_PERM=$n_perm"

Rscript "$POSTDOC_ROOT/OCM_heart/allelic_ratio/13_pseudobulk_celltype_ar.R" || exit 1

out="$results_root/pseudobulk_celltype"
tag="_term${terminal_mb%.0}Mb"
echo
echo "figures:  $PWD/$out/pseudobulk_celltype_chromosome_ends${tag}.pdf"
echo "coverage: $PWD/$out/pseudobulk_end_contrast_coverage${tag}.txt"
echo "  chrX has ~12 genes within 1 Mb of an end and ~60 within 5 Mb, so a 5 Mb"
echo "  window declines most chrX pseudobulks while the autosomes are fine at it."
echo "  Run the X wider - outputs are tagged by window, so nothing is overwritten:"
echo "    sbatch slurm/pseudobulk_celltype_R.slurm analysis $results_root 10"
