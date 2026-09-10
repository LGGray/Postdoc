#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J pseudobulk_celltype_R
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=48G
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
# cm4, not serial_std, even though both stages are single-threaded and small:
# an LRZ --dependency does not cross clusters, and the BAM build and the
# Allelome.PRO2 fan-out are both on cm4. Putting these on serial would break
# the chain, which is the whole point of the script.
#
# 17 cpus for a SINGLE-THREADED job is not a mistake and must not be trimmed.
# cm4_tiny's documented CPU range is 17-112 physical cores -- 17 is a QOS
# FLOOR, not a ceiling, and sbatch rejects anything below it outright with
# "QOSMinCpuNotSatisfied" rather than rounding up. The partition table is at
# https://doku.lrz.de/job-processing-on-the-linux-cluster-10745970.html
# (Do not take a core count from a neighbouring script here: 4 was copied from
# index_igv_bams.slurm, which is itself below the floor and cannot ever have
# been submitted successfully.)
#
# The alternative -- serial_std, whose 1-16 range actually fits a serial job --
# costs the dependency chain, so the trade is 16 idle cores for a few minutes
# against having to babysit four submissions.
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
echo
echo "figures:  $PWD/$out/pseudobulk_celltype_chromosome_ends.pdf"
echo "coverage: $PWD/$out/pseudobulk_end_contrast_coverage.txt"
echo "  If chrX is mostly 'tested FALSE' there, resubmit just this stage with a"
echo "  wider window:  sbatch slurm/pseudobulk_celltype_R.slurm analysis $results_root 10"
