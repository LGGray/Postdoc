#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_signac
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=12:00:00

# Signac: common peak set, cross-sample joint WNN, coverage tracks.
#
# Heavier than multiome_cluster_umap.slurm because of FeatureMatrix, which
# streams both fragment files (723M and 1.1G gzipped) and requantifies every
# QC-passing nucleus against ~100k unified peaks. That is the long pole; expect
# a few hours rather than the ~1h of script 02.
#
# 100G is inside the documented serial_std ceiling (1000 GiB/node - see
# https://doku.lrz.de/job-processing-on-the-linux-cluster-10745970.html) and
# matches the largest memory request already proven in this repo.
#
# seurat_env, which now also carries Signac.

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate seurat_env

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc || exit 1
Rscript multiome/03_signac_joint.R
