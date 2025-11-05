#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH -J RNAseq
#SBATCH -D  /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
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
#SBATCH --array=1-24

module load slurm_setup

#######################################################
# This file is used to preprocess and map RNAseq data #
#######################################################

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

# Creating STAR index - Read length is 50bp so sjdbOverhang is set to 49 ##
STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir mm10_STAR_index \
--genomeFastaFiles /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/RNA_reference_files/genome_ref_mm10/GENCODE_M25GRCm38.p6_201911/GRCm38.primary_assembly.genome.fa \
--sjdbGTFfile /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/RNA_reference_files/genome_ref_mm10/GENCODE_M25GRCm38.p6_201911/gencode.vM25.primary_assembly.annotation.gtf \
--sjdbOverhang 100

## Running STAR alignment - parallel job for the 24 samples ##

# fastq_dir=cardiac_RNAseq/fastq
# RESULTS=cardiac_RNAseq/GRCm39_mapped
# SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" cardiac_RNAseq/fastq_files.txt)


# ## STAR mapping for SNPsplit ##
# STAR \
# --genomeDir  GRCm39/STAR_index \
# --readFilesIn $fastq_dir/${SAMPLE}_1.fastq $fastq_dir/${SAMPLE}_2.fastq \
# --runThreadN 9 \
# --outStd Log \
# --outSAMtype BAM SortedByCoordinate \
# --alignIntronMax 100000 \
# --outFilterIntronMotifs RemoveNoncanonical \
# --outFileNamePrefix $RESULTS/${SAMPLE}"_" \
# --outFilterMultimapNmax 1

# ## Indexing BAM file ##
# samtools index $RESULTS/${SAMPLE}_Aligned.sortedByCoord.out.bam

# ## Couting reads with HTSeq ##
# STRAND=reverse
# htseq-count -s $STRAND -r pos -f bam $RESULTS/$SAMPLE"_Aligned.sortedByCoord.out.bam" GRCm39/gencode.vM37.primary_assembly.annotation.gtf > $RESULTS"/HTseq_counts/"${SAMPLE}"_stranded_"${STRAND}".count"