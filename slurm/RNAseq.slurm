#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH -J RNAseq
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
#SBATCH --cpus-per-task=9
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --array=1-24

module load slurm_setup

#######################################################
# This file is used to preprocess and map RNAseq data #
#######################################################

## Load conda environment ##
conda activate RNAseq

## Creating STAR index - interactive job - salloc -M inter -p cm4_inter -n 1 -c 8 -t 01:00:00 *COMPLETE* ##
# STAR \
# --runThreadN 8 \
# --runMode genomeGenerate \
# --genomeDir genomeDir \
# --genomeFastaFiles genomeDir/GRCm38.primary_assembly.genome.fa \
# --sjdbGTFfile genomeDir/gencode.vM23.annotation.gtf \
# --sjdbOverhang 49

## Running STAR alignment - parallel job for the 24 samples ##

fastq_dir=cardiac_RNAseq
RESULTS=cardiac_RNAseq
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $fastq_dir/samplelist.txt)

STAR \
--genomeDir  genomeDir \
--readFilesIn $fastq_dir/${SAMPLE}_1.fastq $fastq_dir/${SAMPLE}_2.fastq --runThreadN 9 --outStd Log --outSAMtype BAM SortedByCoordinate \
--alignIntronMax 100000 --outFilterIntronMotifs RemoveNoncanonical \
--outFileNamePrefix $RESULTS"/STAR_alignment/"${SAMPLE}"_"  --outFilterMultimapNmax 1

## Indexing BAM file ##
samtools index $RESULTS/STAR_alignment/${SAMPLE}_Aligned.sortedByCoord.out.bam

## Couting reads with HTSeq ##
STRAND=reverse
htseq-count -s $STRAND -r pos -f bam $RESULTS"/STAR_alignment/"$SAMPLE"_Aligned.sortedByCoord.out.bam" genomeDir/gencode.vM23.annotation.gtf > $RESULTS"/HTseq_counts/"${SAMPLE}"_stranded_"${STRAND}".count"