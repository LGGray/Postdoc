#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH -J bowtie
#SBATCH -D  /dss/dssfs02/lwp-dss-0001/pr23fa/pr23fa-dss-0000/go93qiw2
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
#SBATCH --time=1:00:00

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs02/lwp-dss-0001/pr23fa/pr23fa-dss-0000/go93qiw2/miRNA/cardiac_miRNA

# Run bowtie-build
module load bowtie/1.3.1  # Adjust version if needed

bowtie-build -p ${SLURM_CPUS_PER_TASK} ../../genomeDir/GRCm38.primary_assembly.genome.fa GRCm38.primary_assembly.genome