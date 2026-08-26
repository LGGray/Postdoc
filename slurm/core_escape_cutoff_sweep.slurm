#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J core_escape_cutoff_sweep
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

# Cutoff sensitivity sweep for the CORE ESCAPE BLOCK
# (OCM_heart/allelic_ratio/07_core_escape_cutoff_sweep.R).
#
# Needs 04_core_escape.R to have written the raw ratio table already. Cutoffs
# are an order of magnitude lower than the whole-chrX sweep because the block's
# median depth is ~3 reads and its deepest cell has 26.
# Memory is dominated by readRDS on the Seurat object; the sweep itself works
# on a plain data frame of ~1 row per cell.

module load slurm_setup

# --export=NONE means CUTOFFS=... set on the sbatch line does not reach the job.
# Pass the list as the first positional argument instead:
#
#   sbatch core_escape_cutoff_sweep.slurm                 # default sweep
#   sbatch core_escape_cutoff_sweep.slurm 1,3,5,8         # custom cutoffs
#
# Default (set in the R script) is 1,2,3,4,5,8,10.
CUTOFFS="${1:-${CUTOFFS:-}}"

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, not the RNAseq env the rest of this directory uses: the
# allelic_ratio scripts are Seurat-based and need that environment's package
# set. Do not "correct" this to RNAseq for consistency with the neighbours.
conda activate seurat_env

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
cd $BASE/OCM

RATIOS=Allelic_ratio_results/core_escape_block_new_allelic_ratio_table.txt
if [ ! -f "$RATIOS" ]; then
  echo "Missing $RATIOS - run 04_core_escape.R first (it writes this table)." >&2
  exit 1
fi

CUTOFFS="$CUTOFFS" Rscript $BASE/Postdoc/OCM_heart/allelic_ratio/07_core_escape_cutoff_sweep.R
