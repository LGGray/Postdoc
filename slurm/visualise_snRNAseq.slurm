#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J snRNAseq_IGV
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1   
#SBATCH --cpus-per-task=18
#SBATCH --mem=50G
#SBATCH --time=8:00:00
#SBATCH --array=1-10


module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq

celltypes=(Adipo EC Epicard FB Glial LEC Macro SMC atrialCM ventCM)

# index celltypes
ct=${celltypes[$((SLURM_ARRAY_TASK_ID-1))]}

## Adult mice ##
# sort, index
samtools sort -@ 8 -o 9w/sinto/$ct.sorted.bam 9w/sinto/$ct.bam
samtools index 9w/sinto/$ct.sorted.bam

# make WIG from merged BAM
bam2wig.py -i 9w/sinto/$ct.sorted.bam \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes \
  -o 9w/sinto/$ct.adult --skip-multi-hits

## Aged mice ##
# sort, index
samtools sort -@ 8 -o 78w/sinto/$ct.sorted.bam 78w/sinto/$ct.bam
samtools index 78w/sinto/$ct.sorted.bam

# make WIG from merged BAM
bam2wig.py -i 78w/sinto/$ct.sorted.bam \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes \
  -o 78w/sinto/$ct.aged --skip-multi-hits
