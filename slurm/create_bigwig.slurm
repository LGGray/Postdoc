#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J Allelome_LINK
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

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_bodymap

# ##############
# # Adult heart #
# ##############
# # Collect sample names into an array (0-based)
# readarray -t samples < <(grep 'He_9w' adult_metadata.csv | cut -d',' -f2)

# # merge forward
# samtools merge -@ 8 He_9w/pooled_fwd.bam \
#   "adult/Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_fwd.bam" \
#   "adult/Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_fwd.bam" \
#   "adult/Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_fwd.bam"

# samtools sort -o He_9w/pooled_fwd.sorted.bam -@ 8 He_9w/pooled_fwd.bam
# samtools index -@ 8 He_9w/pooled_fwd.sorted.bam

# # star_index=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/STAR2.7.11b_genome_index/genome_dir_gencodeM36_mm39.202410
# # paste $star_index/chrName.txt $star_index/chrLength.txt > ../GRCm39/mm39.chrom.sizes

# bam2wig.py -i He_9w/pooled_fwd.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_9w/He_9w_fwd --skip-multi-hits

# # merge reverse
# samtools merge -@ 8 -f He_9w/pooled_rev.bam \
#   "adult/Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_rev.bam" \
#   "adult/Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_rev.bam" \
#   "adult/Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_rev.bam"

# samtools sort -o He_9w/pooled_rev.sorted.bam -@ 8 He_9w/pooled_rev.bam
# samtools index -@ 8 He_9w/pooled_rev.sorted.bam

# bam2wig.py -i He_9w/pooled_rev.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_9w/He_9w_rev --skip-multi-hits

# # invert counts with deepTools
# bamCoverage \
#   -b He_9w/pooled_rev.sorted.bam \
#   -o He_9w/pooled_rev_inverted.bw \
#   --scaleFactor -1 \
#   --numberOfProcessors 8

# ##############
# # Aged heart #
# ##############
# # Collect sample names into an array (0-based)
# readarray -t samples < <(grep 'He_78w' aged_metadata.csv | cut -d',' -f2)

# # merge forward
# samtools merge -@ 8 He_78w/pooled_fwd.bam \
#   "aged/Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_fwd.bam" \
#   "aged/Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_fwd.bam" \
#   "aged/Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_fwd.bam"

# samtools sort -o He_78w/pooled_fwd.sorted.bam -@ 8 He_78w/pooled_fwd.bam
# samtools index -@ 8 He_78w/pooled_fwd.sorted.bam
# bam2wig.py -i He_78w/pooled_fwd.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_78w/He_78w_fwd --skip-multi-hits

# # merge reverse
# samtools merge -@ 8 He_78w/pooled_rev.bam \
#   "aged/Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_rev.bam" \
#   "aged/Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_rev.bam" \
#   "aged/Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_rev.bam"
# samtools sort -o He_78w/pooled_rev.sorted.bam -@ 8 He_78w/pooled_rev.bam
# samtools index -@ 8 He_78w/pooled_rev.sorted.bam

# bam2wig.py -i He_78w/pooled_rev.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_78w/He_78w_rev --skip-multi-hits

# # invert counts with deepTools
# bamCoverage \
#   -b He_78w/pooled_rev.sorted.bam \
#   -o He_78w/pooled_rev_inverted.bw \
#   --scaleFactor -1 \
#   --numberOfProcessors 8


# #############
# # TAC heart #
# #############

# # Collect sample names into an array (0-based)
# readarray -t samples < <(grep 'He_TAC28d_RNA_XistBxC_-+_XX' metadata.csv | cut -d',' -f2)

# # merge forward
# samtools merge -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.bam \
#   "Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_fwd.bam" \
#   "Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_fwd.bam" \
#   "Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_fwd.bam"

# samtools sort -o He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.sorted.bam -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.bam
# samtools index -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.sorted.bam

# bam2wig.py -i He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_TAC28d_RNA_XistBxC_-+_XX/He_TAC28d_RNA_XistBxC_-+_XX_fwd --skip-multi-hits

# # merge reverse
# samtools merge -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.bam \
#   "Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_rev.bam" \
#   "Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_rev.bam" \
#   "Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_rev.bam"

# samtools sort -o He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.sorted.bam -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.bam
# samtools index -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.sorted.bam

# bam2wig.py -i He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_TAC28d_RNA_XistBxC_-+_XX/He_TAC28d_RNA_XistBxC_-+_XX_rev --skip-multi-hits

# bamCoverage \
#   -b He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.sorted.bam \
#   -o He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev_inverted.bw \
#   --scaleFactor -1 \
#   --numberOfProcessors 8


##############
# Adult FACS #
##############
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/cardiac_RNAseq

# Endothelial cells #
# Collect sample names into an array (0-based)
readarray -t samples < <(grep 'female,Endothelial cells' SraRunTable.csv | cut -d',' -f1)

# merge 
samtools merge -@ 8 female_Endothelial_cells/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"

samtools sort -o female_Endothelial_cells/pooled.sorted.bam -@ 8 female_Endothelial_cells/pooled.bam
samtools index -@ 8 female_Endothelial_cells/pooled.sorted.bam

# star_index=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/STAR2.7.11b_genome_index/genome_dir_gencodeM36_mm39.202410
# paste $star_index/chrName.txt $star_index/chrLength.txt > ../GRCm39/mm39.chrom.sizes

bam2wig.py -i female_Endothelial_cells/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o female_Endothelial_cells/female_Endothelial_cells --skip-multi-hits

# Cardiomyocytes #
# Collect sample names into an array (0-based)
readarray -t samples < <(grep 'female,Cardiomyocytes' SraRunTable.csv | cut -d',' -f1)

# merge
samtools merge -@ 8 female_Cardiomyocytes/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o female_Cardiomyocytes/pooled.sorted.bam -@ 8 female_Cardiomyocytes/pooled.bam
samtools index -@ 8 female_Cardiomyocytes/pooled.sorted.bam

bam2wig.py -i female_Cardiomyocytes/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o female_Cardiomyocytes/female_Cardiomyocytes --skip-multi-hits

# Fibroblasts #
# Collect sample names into an array (0-based)
readarray -t samples < <(grep 'female,Fibroblasts' SraRunTable.csv | cut -d',' -f1)

# merge
samtools merge -@ 8 female_Fibroblasts/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o female_Fibroblasts/pooled.sorted.bam -@ 8 female_Fibroblasts/pooled.bam
samtools index -@ 8 female_Fibroblasts/pooled.sorted.bam

bam2wig.py -i female_Fibroblasts/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o female_Fibroblasts/female_Fibroblasts --skip-multi-hits

# Macrophages #
# Collect sample names into an array (0-based)
readarray -t samples < <(grep 'female,Macrophages' SraRunTable.csv | cut -d',' -f1)

# merge
samtools merge -@ 8 female_Macrophages/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o female_Macrophages/pooled.sorted.bam -@ 8 female_Macrophages/pooled.bam
samtools index -@ 8 female_Macrophages/pooled.sorted.bam

bam2wig.py -i female_Macrophages/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o female_Macrophages/female_Macrophages --skip-multi-hits

#############
# Aged FACS #
#############
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/aged_cardiac_RNAseq

# Endothelial cells #
readarray -t samples < <(grep 'EC_78w_RNA_XBxC_-+_XX' Project_699_lims.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 female_EC/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o female_EC/pooled.sorted.bam -@ 8 female_EC/pooled.bam
samtools index -@ 8 female_EC/pooled.sorted.bam

bam2wig.py -i female_EC/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o female_EC/female_EC --skip-multi-hits

# Cardiomyocytes #
readarray -t samples < <(grep 'CM_78w_RNA_XBxC_-+_XX' Project_699_lims.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 female_CM/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o female_CM/pooled.sorted.bam -@ 8 female_CM/pooled.bam
samtools index -@ 8 female_CM/pooled.sorted.bam

bam2wig.py -i female_CM/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o female_CM/female_CM --skip-multi-hits

# Fibroblasts #
readarray -t samples < <(grep 'FB_78w_RNA_XBxC_-+_XX' Project_699_lims.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 female_FB/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o female_FB/pooled.sorted.bam -@ 8 female_FB/pooled.bam
samtools index -@ 8 female_FB/pooled.sorted.bam

bam2wig.py -i female_FB/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o female_FB/female_FB --skip-multi-hits

# Macrophages #
readarray -t samples < <(grep 'MP_78w_RNA_XBxC_-+_XX' Project_699_lims.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 female_MP/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o female_MP/pooled.sorted.bam -@ 8 female_MP/pooled.bam
samtools index -@ 8 female_MP/pooled.sorted.bam

bam2wig.py -i female_MP/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o female_MP/female_MP --skip-multi-hits

#############
# Sham FACS #
#############
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/TAC_cardiac_RNAseq

# Endothelial cells #
readarray -t samples < <(grep 'EC_Sham_d28_RNA_XBxC' Project_1050_lims_fixed.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 EC_Sham/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o EC_Sham/pooled.sorted.bam -@ 8 EC_Sham/pooled.bam
samtools index -@ 8 EC_Sham/pooled.sorted.bam

bam2wig.py -i EC_Sham/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o EC_Sham/EC_Sham --skip-multi-hits

# Cardiomyocytes #
readarray -t samples < <(grep 'CM_Sham_d28_RNA_XBxC' Project_1050_lims_fixed.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 CM_Sham/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o CM_Sham/pooled.sorted.bam -@ 8 CM_Sham/pooled.bam
samtools index -@ 8 CM_Sham/pooled.sorted.bam

bam2wig.py -i CM_Sham/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o CM_Sham/CM_Sham --skip-multi-hits

# Fibroblasts #
readarray -t samples < <(grep 'FB_Sham_d28_RNA_XBxC' Project_1050_lims_fixed.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 FB_Sham/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o FB_Sham/pooled.sorted.bam -@ 8 FB_Sham/pooled.bam
samtools index -@ 8 FB_Sham/pooled.sorted.bam

bam2wig.py -i FB_Sham/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o FB_Sham/FB_Sham --skip-multi-hits

# Macrophages #
readarray -t samples < <(grep 'MP_Sham_d28_RNA_XBxC' Project_1050_lims_fixed.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 MP_Sham/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o MP_Sham/pooled.sorted.bam -@ 8 MP_Sham/pooled.bam
samtools index -@ 8 MP_Sham/pooled.sorted.bam

bam2wig.py -i MP_Sham/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o MP_Sham/MP_Sham --skip-multi-hits

############
# TAC FACS #
############

# Endothelial cells #
readarray -t samples < <(grep 'EC_TAC_d28_RNA_XBxC' Project_1050_lims_fixed.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 EC_TAC/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o EC_TAC/pooled.sorted.bam -@ 8 EC_TAC/pooled.bam
samtools index -@ 8 EC_TAC/pooled.sorted.bam

bam2wig.py -i EC_TAC/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o EC_TAC/EC_TAC --skip-multi-hits

# Cardiomyocytes #
readarray -t samples < <(grep 'CM_TAC_d28_RNA_XBxC' Project_1050_lims_fixed.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 CM_TAC/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o CM_TAC/pooled.sorted.bam -@ 8 CM_TAC/pooled.bam
samtools index -@ 8 CM_TAC/pooled.sorted.bam

bam2wig.py -i CM_TAC/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o CM_TAC/CM_TAC --skip-multi-hits

# Fibroblasts #
readarray -t samples < <(grep 'FB_TAC_d28_RNA_XBxC' Project_1050_lims_fixed.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 FB_TAC/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o FB_TAC/pooled.sorted.bam -@ 8 FB_TAC/pooled.bam
samtools index -@ 8 FB_TAC/pooled.sorted.bam

bam2wig.py -i FB_TAC/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o FB_TAC/FB_TAC --skip-multi-hits

# Macrophages #
readarray -t samples < <(grep 'MP_TAC_d28_RNA_XBxC' Project_1050_lims_fixed.csv | cut -d',' -f2)

# Merge
samtools merge -@ 8 MP_TAC/pooled.bam \
  "${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out.bam" \
  "${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out.bam" \
  "${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out.bam"
samtools sort -o MP_TAC/pooled.sorted.bam -@ 8 MP_TAC/pooled.bam
samtools index -@ 8 MP_TAC/pooled.sorted.bam
bam2wig.py -i MP_TAC/pooled.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o MP_TAC/MP_TAC --skip-multi-hits

