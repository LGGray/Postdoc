#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_gene_level_escape
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

# Per-gene chrX escape from the gene-level pseudobulk tree. The analysis half
# of slurm/multiome_allelome_genelevel.slurm, which ran the Allelome.PRO2 side
# on 2026-09-11 and has never had its output read.
#
#   sbatch slurm/multiome_gene_level_escape.slurm
#
# ---------------------------------------------------------------------------
# THIS READS DSS DIRECTLY - NOTHING TO COPY
# ---------------------------------------------------------------------------
# Unlike the ATAC job, multiome_allelome_genelevel.slurm already wrote its
# locus tables to adult_aged_multiome/allelome_gene_level/<sample>/, together
# with gene_level_snp_count.bed. So the inputs are on DSS and nothing here
# depends on scratch surviving.
#
# 2h and 48G are generous for 16 locus tables. The time is for the distal
# enrichment block at the end, which sources
# OCM_heart/ap2_pseudobulk_common.R and runs its circular permutation - that is
# the only part of this script that is not near-instant.
#
# ---------------------------------------------------------------------------
# EXPECT THE PRIOR-ESCAPE WARNING TO BE SILENT NOW, AND CHECK IF IT IS NOT
# ---------------------------------------------------------------------------
# 06 warns when the pooled gene-level escape lands more than a factor of two
# from PRIOR_ESCAPE, which is ANALYSIS_PLAN.md step 4's gate. That constant was
# 0.127 and is now 0.029, following the spatial retraction in
# spatial/NEXT_ANALYSIS.md:386-411 - the old figure was chromosome-wide and
# mostly two artefacts. Against the revised prior the multiome lands at roughly
# 0.029 (9w) and 0.040 (78w), so the gate should pass quietly. A warning in the
# log means something moved, not that the old number came back.
#
# Overrides, if they are wanted, go HERE and not in the submitting shell:
# --export=NONE means nothing from that shell arrives.
#
# export MIN_READS=20     # SNP-overlapping reads per gene
# export MIN_SNP=1        # informative SNPs per gene
# ---------------------------------------------------------------------------

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate seurat_env
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc || exit 1
Rscript multiome/06_gene_level_escape.R
