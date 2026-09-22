#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_celltype_violin
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=0:30:00

# Per-cell-type whole-chrX allelic ratio violins, RNA and ATAC as two panels.
#
#   sbatch slurm/multiome_celltype_violin.slurm
#
# ---------------------------------------------------------------------------
# READS CACHES, NOT TREES - WHICH IS WHY THIS IS A 30-MINUTE JOB
# ---------------------------------------------------------------------------
# 09_celltype_ar_violin_rna_atac.R consumes the consolidated per-nucleus TSVs
# that 05_allelome_plots.R and 08_atac_allelome.R already wrote:
#
#   figures_allelome/allelome_pernucleus.tsv        (RNA)
#   figures_atac/allelome_atac_pernucleus.tsv       (ATAC)
#   figures_signac/nucleus_celltype_joint.csv       (cell type labels)
#
# So it never stats the 4,561 Allelome.PRO2 output directories and runs in
# seconds. Both consolidation jobs must have run first; if either cache holds
# only one sample the script STOPS rather than drawing a one-violin figure, and
# the fix is to re-run that job with REBUILD=1.
#
# ---------------------------------------------------------------------------
# WHAT THE DEFAULTS ARE, AND WHEN TO CHANGE THEM
# ---------------------------------------------------------------------------
# Overrides go HERE, not in the submitting shell: --export=NONE means nothing
# from that shell arrives.
#
# export PAIRED=0      # draw every nucleus passing its own gate, instead of
#                      # only the ~2.4k clearing BOTH. The default keeps the
#                      # two panels on one population of cells so the
#                      # comparison between them is paired; PAIRED=0 gains
#                      # ~1.3k RNA and ~440 ATAC nuclei and loses that.
# export MIN_RNA=30    # informative chrX molecules per nucleus
# export MIN_ATAC=50   # informative chrX fragments per nucleus - higher
#                      # because ATAC is ~4x deeper per nucleus here
# export MIN_CELLS=20  # nuclei per cell type x sample, on the WORSE modality
# export MONO_AR=0.90  # OCM_heart/allelic_ratio/00_functions.R; drawn on the
#                      # RNA panel only, since a monoallelic bound is a claim
#                      # about transcription and says nothing about chromatin
#
# Writes celltype_ar_violin_rna_atac.pdf (3 pages) plus two TSVs into
# figures_allelome/, beside 07's single-modality version.
# ---------------------------------------------------------------------------

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate seurat_env
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc || exit 1
Rscript multiome/09_celltype_ar_violin_rna_atac.R
