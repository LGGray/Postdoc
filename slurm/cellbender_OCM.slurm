#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J cellbender_OCM
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
#SBATCH --cpus-per-task=112
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --array=1-4

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate cellbender

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Hearts_OCM/outs/per_sample_outs

## Map array task IDs to input/output paths. Edit or extend as needed.
if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 1 ]; then
	INPUT="He_9w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.h5"
	OUTPUT="He_9w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix_cellbender.h5"
elif [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 2 ]; then
	INPUT="He_78w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.h5"
	OUTPUT="He_78w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix_cellbender.h5"
elif [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 3 ]; then
	INPUT="He_Sham_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.h5"
	OUTPUT="He_Sham_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix_cellbender.h5"
elif [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 4 ]; then
	INPUT="He_TAC_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.h5"
	OUTPUT="He_TAC_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix_cellbender.h5"
else
	echo "Unknown SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}. Exiting." >&2
	exit 1
fi


# Run cellbender with the chosen paths
cellbender remove-background \
	--input "${INPUT}" \
	--output "${OUTPUT}" \
	--expected-cells 6000 \
	--estimator-multiple-cpu

ptrepack --complevel 5 Hearts_OCM/outs/per_sample_outs/He_9w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered.h5:/matrix \
Hearts_OCM/outs/per_sample_outs/He_9w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered_seurat.h5:/matrix

ptrepack --complevel 5 Hearts_OCM/outs/per_sample_outs/He_78w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered.h5:/matrix \
Hearts_OCM/outs/per_sample_outs/He_78w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered_seurat.h5:/matrix

ptrepack --complevel 5 Hearts_OCM/outs/per_sample_outs/He_Sham_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered.h5:/matrix \
Hearts_OCM/outs/per_sample_outs/He_Sham_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered_seurat.h5:/matrix

ptrepack --complevel 5 Hearts_OCM/outs/per_sample_outs/He_TAC_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered.h5:/matrix \
Hearts_OCM/outs/per_sample_outs/He_TAC_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered_seurat.h5:/matrix