#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J sc_Allelome.PRO2.ShamTAC
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

# ### 9w
# # Get list of completed barcodes from existing locus_table paths
# ls Allelome.PRO2/9w/*/locus_table.txt | \
#   sed 's|.*/2026_06_12_9w_\(.*\)_chr_annotation.*|\1|' > Allelome.PRO2/9w/completed.txt

# # Extract barcode from task-list (the BAM filename without .bam)
# grep -oP '9w_sc_sinto/\K[^.]+(?=\.bam)' Allelome.PRO2/9w/task-list_noXist.txt > /tmp/all.txt

# # Extract barcode from completed locus_tables
# ls Allelome.PRO2/9w/*/locus_table.txt | grep -oP '9w_\K[^_chr]+(?=_chr)' > /tmp/completed.txt

# # Check counts look right
# wc -l /tmp/all.txt /tmp/completed.txt

# # Add 9w_ prefix to completed list
# sed -i 's/^/9w_/' /tmp/completed.txt

# # Verify they match now
# head -3 /tmp/completed.txt

# # Check counts
# wc -l /tmp/all.txt /tmp/completed.txt

# # Generate remaining task-list
# grep -Fvf /tmp/completed.txt Allelome.PRO2/9w/task-list_noXist.txt \
#   > Allelome.PRO2/9w/task-list_remaining.txt

# wc -l Allelome.PRO2/9w/task-list_remaining.txt

parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2/9w/task-list_remaining.txt

# # ### 78w
# # Extract barcode from task-list
# grep -oP '78w_sc_sinto/\K[^.]+(?=\.bam)' Allelome.PRO2/78w/task-list_noXist.txt > /tmp/all_78w.txt

# # Extract barcode from completed locus_tables
# ls Allelome.PRO2/78w/*/locus_table.txt | grep -oP '78w_\K[^_chr]+(?=_chr)' > /tmp/completed_78w.txt

# # Add prefix
# sed -i 's/^/78w_/' /tmp/completed_78w.txt

# # Check counts
# wc -l /tmp/all_78w.txt /tmp/completed_78w.txt

# # Generate remaining task-list
# grep -Fvf /tmp/completed_78w.txt Allelome.PRO2/78w/task-list_noXist.txt \
#   > Allelome.PRO2/78w/task-list_remaining.txt

# wc -l Allelome.PRO2/78w/task-list_remaining.txt

parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2/78w/task-list_remaining.txt

# ######------ Sham ------######
pipeline="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh"
input="sinto/Sham/bamfiles.txt"
annotation="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/chr_annotation_mm39.bed"
snps="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.clean.bed"
output_dir="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2/Sham"
min_read=1
total_read=1

while IFS= read -r line; do
  echo "$pipeline -i sinto/Sham/$line -a $annotation -s $snps -r $min_read -t $total_read -o $output_dir"
  done < "$input" > Allelome.PRO2/Sham/task-list_noXist.txt

parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2/Sham/task-list_noXist.txt


# ######------ TAC ------######
pipeline="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh"
input="sinto/TAC/bamfiles.txt"
annotation="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/chr_annotation_mm39.bed"
snps="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.clean.bed"
output_dir="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2/TAC"
min_read=1
total_read=1

while IFS= read -r line; do
  echo "$pipeline -i sinto/TAC/$line -a $annotation -s $snps -r $min_read -t $total_read -o $output_dir"
  done < "$input" > Allelome.PRO2/TAC/task-list_noXist.txt

parallel -j "$nproc" :::: /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM/Allelome.PRO2/TAC/task-list_noXist.txt