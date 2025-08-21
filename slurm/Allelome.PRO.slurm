#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH -J Allelome.PRO
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
#SBATCH --array=1-24

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" cardiac_RNAseq/fastq_files.txt)

bash Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/cardiac_RNAseq/GRCm39_mapped/${SAMPLE}_Aligned.sortedByCoord.out.bam \
  -a /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/cardiac_RNAseq/Allelome.Pro_GRCm39/20250512_RefSeq.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/cardiac_RNAseq/Allelome.Pro_GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_exon_unique_REL2021.bed \
  -o /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/cardiac_RNAseq/Allelome.Pro_GRCm39