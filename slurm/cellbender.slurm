#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J cellbender
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
#SBATCH --cpus-per-task=18
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --array=1-2

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate cellbender

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq

## Map array task IDs to input/output paths. Edit or extend as needed.
if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 1 ]; then
	INPUT="9w/outs/raw_feature_bc_matrix.h5"
	OUTPUT="9w/outs/raw_feature_bc_matrix_cellbender.h5"
elif [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 2 ]; then
	INPUT="78w/outs/raw_feature_bc_matrix.h5"
	OUTPUT="78w/outs/raw_feature_bc_matrix_cellbender.h5"
else
	echo "Unknown SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}. Exiting." >&2
	exit 1
fi

# Run cellbender with the chosen paths
cellbender remove-background \
	--input "${INPUT}" \
	--output "${OUTPUT}"
