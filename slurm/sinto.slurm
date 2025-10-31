#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J sinto
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
#SBATCH --mem=50G
#SBATCH --time=8:00:00
#SBATCH --array=1-2

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate sinto

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq

## Map array task IDs to input/output paths. Edit or extend as needed.
if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 1 ]; then
	input_bam="9w/outs/possorted_genome_bam.bam"
	out_cell_barcodes_annot="9w/cell_index.txt"
    output_dir="9w/sinto"
elif [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 2 ]; then
	input_bam="78w/outs/possorted_genome_bam.bam"
	out_cell_barcodes_annot="78w/cell_index.txt"
    output_dir="78w/sinto"
else
	echo "Unknown SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}. Exiting." >&2
	exit 1
fi

mkdir -p "$output_dir"
# Use the SLURM allocation for CPUs if available, otherwise fall back to 8
nproc=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS="$nproc"
sinto filterbarcodes --bam "$input_bam" --cells "$out_cell_barcodes_annot" --outdir "$output_dir" --nproc "$nproc"

if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 1 ]; then
	input_bam="9w/outs/possorted_genome_bam.bam"
	out_cell_barcodes_annot="pseudobulk/adult_cell_index.txt"
    output_dir="pseudobulk/"
elif [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 2 ]; then
	input_bam="78w/outs/possorted_genome_bam.bam"
	out_cell_barcodes_annot="pseudobulk/aged_cell_index.txt"
	output_dir="pseudobulk/"
else
	echo "Unknown SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}. Exiting." >&2
	exit 1
fi

nproc=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS="$nproc"
sinto filterbarcodes --bam "$input_bam" --cells "$out_cell_barcodes_annot" --outdir "$output_dir" --nproc "$nproc"