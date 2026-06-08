#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J sc_Allelome.PRO2.snRNAseq
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --cpus-per-task=56
#SBATCH --mem=240G
#SBATCH --ntasks=1
#SBATCH --qos=cm4_tiny
#SBATCH --time=12:00:00
#SBATCH --export=NONE

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM

# ######------ make task-list ------######
# #generate the task-db file
# pipeline="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh"
# input="$SCRATCH/sinto_180476/bamfiles.txt"
# annotation="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/chr_annotation_mm39.bed"
# snps="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.bed"
# output_dir="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2/9w"
# min_read=1
# total_read=1

# while IFS= read -r line; do
#   echo "$pipeline -i $SCRATCH/sinto_180476/$line -a $annotation -s $snps -r $min_read -t $total_read -o $output_dir"
#   done < "$input" > Allelome.PRO2/9w/task-list_noXist.txt

nproc=$((SLURM_CPUS_PER_TASK > 10 ? SLURM_CPUS_PER_TASK - 10 : 1))
module load parallel
parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2/9w/task-list_noXist.txt
