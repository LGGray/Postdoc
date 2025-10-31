#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J Allelome.LINK.snRNAseq
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1   
#SBATCH --cpus-per-task=18
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --array=1-20

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filelist.txt)
AGE=$(echo $SAMPLE | cut -d/ -f1)
celltype=$(echo $SAMPLE | cut -d/ -f3 | sed 's/\.bam$//')

# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i $SAMPLE\
#   -a /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/20250512_RefSeq.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/$AGE/Allelome.LINK

# out_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/${AGE}/Allelome.LINK/2025_09_08_${celltype}_20250512_RefSeq.bed_1/
# mkdir -p ${out_dir}/LINK

# Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
#   -i ${out_dir}/locus_table.txt \
#   -o ${out_dir}/LINK/

# SAMPLE=9w/sinto/Epicard.bam
# AGE=9w
# celltype=Epicard

ref_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39

# Allelome Pro with no predicted genes
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i $SAMPLE\
  -a $ref_dir/annotation_no_predicted_us.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
  -o /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/$AGE/Allelome.LINK

# # Allelome Pro with predicted genes
# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i $SAMPLE\
#   -a $ref_dir/annotation_us.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/$AGE/Allelome.LINK

DATE=$(date +%Y_%m_%d)
out_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/${AGE}/Allelome.LINK/${DATE}_${celltype}_annotation_no_predicted_us.bed_1

Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
  -i ${out_dir}/locus_table.txt \
  -o ${out_dir} \
  --allelic-bias 0.7 \
  -n ${celltype}_no_predicted_70

Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
  -i ${out_dir}/locus_table.txt \
  -o ${out_dir} \
  --allelic-bias 0.65 \
  -n ${celltype}_no_predicted_65

# out_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/${AGE}/Allelome.LINK/${DATE}_${celltype}_annotation_us.bed_1
# mkdir -p ${out_dir}/LINK

# Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
#   -i ${out_dir}/locus_table.txt \
#   -o ${out_dir} \
#   --allelic-bias 0.70 \
#   -n ${celltype}_70

# Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
#   -i ${out_dir}/locus_table.txt \
#   -o ${out_dir} \
#   --allelic-bias 0.65 \
#   -n ${celltype}_65