#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH -J adult_Allelome_lINK
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
#SBATCH --array=1-20

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

# Change to the working directory
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_bodymap

# Get sample name
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" adult_path.txt)
# Filename is SAMPLE split by _ and take the second part
filename=$(echo $SAMPLE | cut -d'_' -f2)

# make a directory for each sample
mkdir -p adult
mkdir -p adult/${SAMPLE}

adult_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/01_aging_TH_for_LG/00_rawdata/01_adult

R1_files=(${adult_dir}/${SAMPLE}/${filename}_S*_L*_R1_001.fastq.gz)
R2_files=(${adult_dir}/${SAMPLE}/${filename}_S*_L*_R2_001.fastq.gz)


# collect and sort lanes to keep L001 before L002
shopt -s nullglob
R1_files=($(printf "%s\n" ${adult_dir}/${SAMPLE}/${filename}_S*_L*_R1_001.fastq.gz | sort))
R2_files=($(printf "%s\n" ${adult_dir}/${SAMPLE}/${filename}_S*_L*_R2_001.fastq.gz | sort))
shopt -u nullglob

# join arrays with commas (one string per mate)
IFS=, R1_join="${R1_files[*]}"; unset IFS
IFS=, R2_join="${R2_files[*]}"; unset IFS

STAR \
  --genomeDir ../GRCm39/STAR_index \
  --readFilesIn "$R1_join" "$R2_join" \
  --readFilesCommand zcat \
  --runThreadN 9 \
  --outStd Log \
  --outSAMtype BAM SortedByCoordinate \
  --alignIntronMax 100000 \
  --outFilterIntronMotifs RemoveNoncanonical \
  --outFileNamePrefix "adult/${SAMPLE}/${filename}_" \
  --outFilterMultimapNmax 1

## Indexing BAM file ##
samtools index adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam

## Couting reads with HTSeq ##
mkdir -p adult/${SAMPLE}/HTseq_counts
STRAND=reverse
htseq-count -s $STRAND -r pos -f bam adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam ../GRCm39/gencode.vM37.primary_assembly.annotation.gtf > adult/${SAMPLE}/HTseq_counts/${filename}_stranded_${STRAND}.count
