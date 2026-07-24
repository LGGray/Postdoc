#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J sc_Allelome.PRO2.chrX_gene.9w
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --cpus-per-task=110
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
pipeline="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh"
input="sinto/9w/bamfiles.txt"
annotation="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/annotation_us_mm39_chrX.bed"
snps="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed"
scratch_dir="$SCRATCH/9w"
final_dir="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2_all_genes/9w"
min_read=1
total_read=1

mkdir -p "$scratch_dir"

while IFS= read -r line; do
  echo "$pipeline -i sinto/9w/$line -a $annotation -s $snps -r $min_read -t $total_read -o $scratch_dir"
  done < "$input" > Allelome.PRO2/9w/task-list_noXist.txt

parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2/9w/task-list_noXist.txt

# move locus_table.txt and BED_files/*.bed from scratch into permanent storage, one destination dir per sample
for sample_dir in "$scratch_dir"/*/; do
  sample_name=$(basename "$sample_dir")
  dest_dir="$final_dir/$sample_name"
  mkdir -p "$dest_dir"
  mv "$sample_dir"locus_table.txt "$dest_dir/"
  mv "$sample_dir"BED_files/*.bed "$dest_dir/"
done
