#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J sc_Allelome.PRO2.core_escape
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --ntasks=1
#SBATCH --qos=cm4_tiny
#SBATCH --time=24:00:00
#SBATCH --export=NONE


module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM

nproc=$((SLURM_CPUS_PER_TASK - 2))
if [ "$nproc" -lt 1 ]; then nproc=1; fi

module load parallel

# ######------ 9w ------######
# pipeline="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh"
# input="sinto/9w/bamfiles.txt"
# annotation="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/annotation_us_mm39_chrX.bed"
# snps="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.clean.bed"
# output_dir="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2_all_genes/9w"
# min_read=1
# total_read=1

# while IFS= read -r line; do
#   echo "$pipeline -i sinto/9w/$line -a $annotation -s $snps -r $min_read -t $total_read -o $output_dir"
#   done < "$input" > Allelome.PRO2_all_genes/9w/task-list_noXist.txt

parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2_all_genes/9w/task-list_noXist.txt

# ######------ 78w ------######
# pipeline="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh"
# input="sinto/78w/bamfiles.txt"
# annotation="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/annotation_us_mm39_chrX.bed"
# snps="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.clean.bed"
# output_dir="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2_all_genes/78w"
# min_read=1
# total_read=1

# while IFS= read -r line; do
#   echo "$pipeline -i sinto/78w/$line -a $annotation -s $snps -r $min_read -t $total_read -o $output_dir"
#   done < "$input" > Allelome.PRO2_all_genes/78w/task-list_noXist.txt

parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2_all_genes/78w/task-list_noXist.txt

# ######------ Sham ------######
# pipeline="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh"
# input="sinto/Sham/bamfiles.txt"
# annotation="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/annotation_us_mm39_chrX.bed"
# snps="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.clean.bed"
# output_dir="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2_all_genes/Sham"
# min_read=1
# total_read=1

# while IFS= read -r line; do
#   echo "$pipeline -i sinto/Sham/$line -a $annotation -s $snps -r $min_read -t $total_read -o $output_dir"
#   done < "$input" > Allelome.PRO2_all_genes/Sham/task-list_noXist.txt

parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2_all_genes/Sham/task-list_noXist.txt

# ######------ TAC ------######
# pipeline="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh"
# input="sinto/TAC/bamfiles.txt"
# annotation="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/annotation_us_mm39_chrX.bed"
# snps="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.clean.bed"
# output_dir="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2_all_genes/TAC"
# min_read=1
# total_read=1

# while IFS= read -r line; do
#   echo "$pipeline -i sinto/TAC/$line -a $annotation -s $snps -r $min_read -t $total_read -o $output_dir"
#   done < "$input" > Allelome.PRO2_all_genes/TAC/task-list_noXist.txt

parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2_all_genes/TAC/task-list_noXist.txt