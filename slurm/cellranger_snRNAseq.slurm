#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH -J download_SRA
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_std
#SBATCH --qos=cm4_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=480G
#SBATCH --time=8:00:00
#SBATCH --array=1-2

export PATH=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/cellranger-9.0.1:$PATH

REF=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/cellranger-9.0.1/refdata-gex-GRCm39-2024-A/

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq

# Aged
cellranger count \
  --id=78w \
  --transcriptome="$REF" \
  --fastqs=Aged \
  --sample=Aged \
  --include-introns=true \
  --create-bam=true \
  --localcores=8 --localmem=64

# # Adult
# cellranger count \
#   --id=adult_aged_heart_snRNAseq/Adult \
#   --transcriptome="$REF" \
#   --fastqs="$RAW/Adult" \
#   --sample=Adult \
#   --include-introns=true \
#   --create-bam=true \
#   --localcores=8 --localmem=64