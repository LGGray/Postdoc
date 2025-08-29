#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH -J download_SRA
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

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
export PATH=cellranger-9.0.1:$PATH

# Corrected sample lists (fixed typo SRR30223486)
adult=(SRR30223485 SRR30223486)
aged=(SRR30223487 SRR30223488)

cellranger=cellranger-9.0.1/
reference=$cellranger/refdata-gex-GRCm39-2024-A

# FASTQ_DIR should point to the directory on the cluster where your
# per-SRR FASTQ subdirectories or files are located. Update as needed.
FASTQ_DIR=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq

# Run cellranger count per SRR. Cell Ranger will pick up paired files
# (e.g. SRR30223485_1.fastq.gz and SRR30223485_2.fastq.gz) and multiple
# lane files with the same SRR prefix if present in $FASTQ_DIR.
for s in "${adult[@]}"; do
  cellranger count \
    --id="Adult_${s}" \
    --transcriptome=$reference \
    --fastqs=$FASTQ_DIR \
    --sample=$s \
    --include-introns \
    --localcores=8 --localmem=64
done

for s in "${aged[@]}"; do
  cellranger count \
    --id="Aged_${s}" \
    --transcriptome=$reference \
    --fastqs=$FASTQ_DIR \
    --sample=$s \
    --include-introns \
    --localcores=8 --localmem=64
done