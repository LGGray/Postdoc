#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_allelome_plots
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
#SBATCH --time=3:00:00

# Consolidate and plot whatever multiome_allelome_run.slurm finished.
#
# Submit with afterANY, not afterok:
#   RUN=$(sbatch --parse slurm/multiome_allelome_run.slurm)
#   sbatch --dependency=afterany:${RUN%%;*} slurm/multiome_allelome_plots.slurm
# or just use submit_allelome_chain.sh. afterany matters: if the run job hits
# the 24h wall clock partway through the per-nucleus stage, we still want the
# per-celltype results and however many nuclei finished to be plotted, rather
# than nothing.

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate seurat_env
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc || exit 1
Rscript multiome/05_allelome_plots.R
