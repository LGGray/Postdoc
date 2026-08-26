#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J lox_calling
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
#SBATCH --time=3:00:00

# Does "AR >= 0.90" actually mean a cell has lost its inactive X?
# (OCM_heart/allelic_ratio/08_lox_calling.R)
#
# Needs 04_core_escape.R for the raw core-escape ratio table, and ideally
# 06_cutoff_sweep.R for the whole-chromosome cross-check (that step is skipped
# with a message if its table is absent). Reads the Seurat object for Xist.
#
# Writes a plain-language summary to
# Allelic_ratio_results/lox_calling/00_SUMMARY.txt -- read that first.
# Memory is dominated by readRDS on the Seurat object; the sweep itself works
# on a plain data frame of ~1 row per cell.

module load slurm_setup


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

Rscript $BASE/Postdoc/OCM_heart/allelic_ratio/08_lox_calling.R
