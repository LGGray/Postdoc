#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J Allelome.LINK.pseudobulk
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
#SBATCH --array=1-2

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq

ref_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39
if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 1 ]; then
    bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
    -i pseudobulk/adult.bam \
    -a $ref_dir/annotation_no_predicted_us.bed \
    -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
    -o /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/pseudobulk
elif [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 2 ]; then
    bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
    -i pseudobulk/aged.bam \
    -a $ref_dir/annotation_no_predicted_us.bed \
    -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
    -o /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/pseudobulk
fi

DATE=$(date +%Y_%m_%d)
if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 1 ]; then
    out_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/pseudobulk/${DATE}_adult_annotation_no_predicted_us.bed_1
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
elif [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 2 ]; then
    out_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/pseudobulk/${DATE}_aged_annotation_no_predicted_us.bed_1
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
fi
