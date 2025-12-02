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

##############
# Adult heart #
##############
# Collect sample names into an array (0-based)
readarray -t samples < <(grep 'He_9w' adult_metadata.csv | cut -d',' -f2)

# merge forward
samtools merge -@ 8 He_9w/pooled_fwd.bam \
  "adult/Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_fwd.bam" \
  "adult/Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_fwd.bam" \
  "adult/Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_fwd.bam"

samtools sort -o He_9w/pooled_fwd.sorted.bam -@ 8 He_9w/pooled_fwd.bam
samtools index -@ 8 He_9w/pooled_fwd.sorted.bam

# star_index=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/STAR2.7.11b_genome_index/genome_dir_gencodeM36_mm39.202410
# paste $star_index/chrName.txt $star_index/chrLength.txt > ../GRCm39/mm39.chrom.sizes

bam2wig.py -i He_9w/pooled_fwd.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_9w/He_9w_fwd --skip-multi-hits

# merge reverse
samtools merge -@ 8 He_9w/pooled_rev.bam \
  "adult/Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_rev.bam" \
  "adult/Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_rev.bam" \
  "adult/Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_rev.bam"
samtools sort -o He_9w/pooled_rev.sorted.bam -@ 8 He_9w/pooled_rev.bam
samtools index -@ 8 He_9w/pooled_rev.sorted.bam

bam2wig.py -i He_9w/pooled_rev.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_9w/He_9w_rev --skip-multi-hits

##############
# Aged heart #
##############
# Collect sample names into an array (0-based)
readarray -t samples < <(grep 'He_78w' aged_metadata.csv | cut -d',' -f2)

# merge forward
samtools merge -@ 8 He_78w/pooled_fwd.bam \
  "aged/Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_fwd.bam" \
  "aged/Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_fwd.bam" \
  "aged/Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_fwd.bam"

samtools sort -o He_78w/pooled_fwd.sorted.bam -@ 8 He_78w/pooled_fwd.bam
samtools index -@ 8 He_78w/pooled_fwd.sorted.bam
bam2wig.py -i He_78w/pooled_fwd.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_78w/He_78w_fwd --skip-multi-hits

# merge reverse
samtools merge -@ 8 He_78w/pooled_rev.bam \
  "aged/Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_rev.bam" \
  "aged/Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_rev.bam" \
  "aged/Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_rev.bam"
samtools sort -o He_78w/pooled_rev.sorted.bam -@ 8 He_78w/pooled_rev.bam
samtools index -@ 8 He_78w/pooled_rev.sorted.bam

bam2wig.py -i He_78w/pooled_rev.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_78w/He_78w_rev --skip-multi-hits




# Collect sample names into an array (0-based)
readarray -t samples < <(grep 'He_TAC28d_RNA_XistBxC_-+_XX' metadata.csv | cut -d',' -f2)

# merge forward
samtools merge -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.bam \
  "Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_fwd.bam" \
  "Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_fwd.bam" \
  "Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_fwd.bam"

samtools sort -o He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.sorted.bam -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.bam
samtools index -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.sorted.bam

bam2wig.py -i He_TAC28d_RNA_XistBxC_-+_XX/pooled_fwd.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_TAC28d_RNA_XistBxC_-+_XX/He_TAC28d_RNA_XistBxC_-+_XX_fwd --skip-multi-hits

# merge reverse
samtools merge -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.bam \
  "Sample_${samples[0]}/${samples[0]}_Aligned.sortedByCoord.out_rev.bam" \
  "Sample_${samples[1]}/${samples[1]}_Aligned.sortedByCoord.out_rev.bam" \
  "Sample_${samples[2]}/${samples[2]}_Aligned.sortedByCoord.out_rev.bam"

samtools sort -o He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.sorted.bam -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.bam
samtools index -@ 8 He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.sorted.bam

bam2wig.py -i He_TAC28d_RNA_XistBxC_-+_XX/pooled_rev.sorted.bam -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/mm39.chrom.sizes -o He_TAC28d_RNA_XistBxC_-+_XX/He_TAC28d_RNA_XistBxC_-+_XX_rev --skip-multi-hits

