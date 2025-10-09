#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J testing_code
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1   
#SBATCH --cpus-per-task=18
#SBATCH --mem=8G
#SBATCH --time=1:00:00

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_bodymap

adult_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/01_aging_TH_for_LG/00_rawdata/01_adult


i=1
SAMPLE=$(sed -n "${i}p" adult_path.txt)
filename=$(echo "$SAMPLE" | cut -d'_' -f2)

# Merge locus_table from forward and reverse strands
Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R adult SAMPLE

# Merge across replicates
#Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_replicates.R

# Allelome LINK
#Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
#  -i adult/${SAMPLE}/locus_table_reps.txt \
#  -o ${out_dir}/LINK/