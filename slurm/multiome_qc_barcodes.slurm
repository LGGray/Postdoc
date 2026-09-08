#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_qc
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00

# Multiome QC and per-nucleus barcode export. Small job: it reads the two
# FILTERED matrices (1758 and 3822 nuclei) plus per_barcode_metrics.csv, not
# the raw matrices or the BAMs.
#
# RNA cutoffs come from ddqcR, which clusters internally and takes MAD-based
# cutoffs within each cluster -- hence the headroom above, and hence the choice:
# cardiomyocytes are the most mitochondria-rich cell type in the heart, so any
# GLOBAL mito cutoff preferentially deletes the cell type this project is
# about. The inherited fixed percent.mt < 5 removed 84% of nuclei on its own.
# ATAC floors are applied separately, since ddqcR is RNA-only.
#
# This is the gate on the per-nucleus allelic work. Its output fixes which
# nuclei everything downstream is scored on, so review
# adult_aged_multiome/qc/ before splitting 5580 BAMs.
#
# seurat_env, NOT the RNAseq env that Allelome.PRO2 needs. The allelic scripts
# in this project are split across the two and the neighbouring slurm files are
# misleading about which is which.

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate seurat_env

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
cd "$BASE/Postdoc" || exit 1

Rscript multiome/01_qc_barcodes.R
