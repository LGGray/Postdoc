#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J allelic_ratio_cutoff_sweep
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=2:00:00

# Total-read cutoff sensitivity sweep for the whole-chrX allelic ratio
# (OCM_heart/allelic_ratio/06_cutoff_sweep.R). Re-makes the AR UMAP and the
# per-celltype violins at a range of total_reads thresholds.
#
# Needs 02_whole_chrX.R to have written whole_chr_allelic_ratios.txt already.
# Memory is dominated by readRDS on the Seurat object; the sweep itself works
# on a plain data frame of ~1 row per cell.

module load slurm_setup

# --export=NONE means CUTOFFS=... set on the sbatch line does not reach the job.
# Pass the list as the first positional argument instead:
#
#   sbatch allelic_ratio_cutoff_sweep.slurm                 # default sweep
#   sbatch allelic_ratio_cutoff_sweep.slurm 1,5,10,20,50    # custom cutoffs
#
# Default (set in the R script) is 1,5,10,20,30,50,100.
CUTOFFS="${1:-${CUTOFFS:-}}"

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
cd $BASE/OCM_heart

RATIOS=Allelic_ratio_results/whole_chr_allelic_ratios.txt
if [ ! -f "$RATIOS" ]; then
  echo "Missing $RATIOS - run 02_whole_chrX.R first (it writes this table)." >&2
  exit 1
fi

CUTOFFS="$CUTOFFS" Rscript $BASE/Postdoc/OCM_heart/allelic_ratio/06_cutoff_sweep.R
