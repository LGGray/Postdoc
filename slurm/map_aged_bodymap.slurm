#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J aged_Allelome_LINK
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
#SBATCH --array=1-29


module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

# Change to the working directory
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_bodymap

aged_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/01_aging_TH_for_LG/00_rawdata/02_aged

run_one() {
  SAMPLE="$1"
  filename=$(echo "$SAMPLE" | cut -d'_' -f2)

  mkdir -p "aged/${SAMPLE}"

  R1_files=(${aged_dir}/${SAMPLE}/${filename}_S*_L*_R1_001.fastq.gz)
  R2_files=(${aged_dir}/${SAMPLE}/${filename}_S*_L*_R2_001.fastq.gz)

  # One STAR per srun; --exclusive ensures each gets its own CPU set
  srun --exclusive -n 1 -c "${SLURM_CPUS_PER_TASK}" \
    STAR \
      --genomeDir ../GRCm39/STAR_index \
      --readFilesIn "$R1_files" "$R2_files" \
      --readFilesCommand zcat \
      --runThreadN "${SLURM_CPUS_PER_TASK}" \
      --outStd Log \
      --outSAMtype BAM SortedByCoordinate \
      --alignIntronMax 100000 \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outFileNamePrefix "aged/${SAMPLE}/${filename}_" \
      --outFilterMultimapNmax 1

  samtools index "aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam"

  STRAND=reverse
  htseq-count -s $STRAND -r pos -f bam \
    "aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam" \
    ../GRCm39/gencode.vM37.primary_assembly.annotation.gtf \
    > "aged/${SAMPLE}/${filename}_stranded_${STRAND}.count"

  echo "Finished processing sample: $i $SAMPLE"
}

i=${SLURM_ARRAY_TASK_ID}
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" aged_path.txt)
filename=$(echo "$SAMPLE" | cut -d'_' -f2)
# run_one "$SAMPLE"

STRAND=reverse
htseq-count -s $STRAND -r pos -f bam \
  "aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam" \
  ../GRCm39/gencode.vM37.primary_assembly.annotation.nochr.gtf \
  > "aged/${SAMPLE}/${filename}_stranded_${STRAND}.count"

# # # Separate BAM files into forward and reverse reads
# strandrule="1+-,1-+,2++,2--"
# perl ../Postdoc/separate_BAM_strand.pl "aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam" $strandrule aged/${SAMPLE}

ref_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39

# Allelome Pro with no predicted genes
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
  -a $ref_dir/annotation_no_predicted_f.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
  -o "aged/${SAMPLE}"

bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_rev.bam \
  -a $ref_dir/annotation_no_predicted_r.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
  -o "aged/${SAMPLE}"

# Allelome Pro with predicted genes
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
  -a $ref_dir/annotation_f.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
  -o "aged/${SAMPLE}"

bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_rev.bam \
  -a $ref_dir/annotation_r.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
  -o "aged/${SAMPLE}"


# Merge locus_table from forward and reverse strands
DATE=$(date +%Y_%m_%d)
Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R aged/$SAMPLE $DATE annotation_no_predicted
Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R aged/$SAMPLE $DATE annotation