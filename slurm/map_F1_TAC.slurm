#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J F1_TAC
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
#SBATCH --array=1-12

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/F1_TAC_Sarah

run_one() {
  SAMPLE="$1"
  filename=$(echo "$SAMPLE" | cut -d'_' -f2)

  R1_files=(${SAMPLE}/${filename}_S*_L*_R1_001.fastq.gz)
  R2_files=(${SAMPLE}/${filename}_S*_L*_R2_001.fastq.gz)

  # One STAR per srun; --exclusive ensures each gets its own CPU set
  srun --exclusive -n 1 -c "${SLURM_CPUS_PER_TASK}" \
    STAR \
      --genomeDir /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/STAR_index \
      --readFilesIn "$R1_files" "$R2_files" \
      --readFilesCommand zcat \
      --runThreadN "${SLURM_CPUS_PER_TASK}" \
      --outStd Log \
      --outSAMtype BAM SortedByCoordinate \
      --alignIntronMax 100000 \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outFileNamePrefix "${SAMPLE}/${filename}_" \
      --outFilterMultimapNmax 1

  samtools index "${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam"

  STRAND=reverse
  htseq-count -s $STRAND -r pos -f bam \
    "${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam" \
    /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/gencode.vM37.primary_assembly.annotation.gtf \
    > "${SAMPLE}/${filename}_stranded_${STRAND}.count"

  echo "Finished processing sample: $i $SAMPLE"
}


SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_names.txt)
filename=$(echo "$SAMPLE" | cut -d'_' -f2)
# run_one "$SAMPLE"

# STRAND=reverse
# htseq-count -s $STRAND -r pos -f bam \
#   "${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam" \
#   /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/gencode.vM37.primary_assembly.annotation.nochr.gtf \
#   > "${SAMPLE}/${filename}_stranded_${STRAND}.count"

# # Separate BAM files into forward and reverse reads
#BAM="${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam"
#ANNOTATION=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/gencode.vM37.primary_assembly.annotation.gtf
#BED=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/gencode.vM37.primary_assembly.annotation.bed

#python /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/infer_experiment.py -i $BAM -r $BED

# strandrule="1+-,1-+,2++,2--"
# perl /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/separate_BAM_strand.pl "${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam" $strandrule ${SAMPLE}

ref_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39

# Allelome Pro with no predicted genes
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i ${SAMPLE}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
  -a $ref_dir/annotation_no_predicted_f.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
  -o "${SAMPLE}"

bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i ${SAMPLE}/${filename}_Aligned.sortedByCoord.out_rev.bam \
  -a $ref_dir/annotation_no_predicted_r.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
  -o "${SAMPLE}"

# Allelome Pro with predicted genes
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i ${SAMPLE}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
  -a $ref_dir/annotation_f.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
  -o "${SAMPLE}"

bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i ${SAMPLE}/${filename}_Aligned.sortedByCoord.out_rev.bam \
  -a $ref_dir/annotation_r.bed \
  -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
  -o "${SAMPLE}"

# Merge locus_table from forward and reverse strands
DATE=$(date +%Y_%m_%d)
Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R $SAMPLE $DATE annotation_no_predicted
Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R $SAMPLE $DATE annotation

# for i in {1..12}; do
#   SAMPLE=$(sed -n "${i}p" sample_names.txt)
#   Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R $SAMPLE $DATE annotation_no_predicted
#   Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R $SAMPLE $DATE annotation
# done

# ##Allelome LINK analysis for predicted and no predicted genes at 65% and 70% ratios
# for i in {1..4}; do
#   SAMPLE=$(sed -n "${i}p" dir_names.txt)
#   rm -r ${SAMPLE}/no_predicted_70
#   Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
#   -i ${SAMPLE}/annotation_no_predicted_locus_table_reps.txt \
#   -o ${SAMPLE} \
#   --allelic-bias 0.7 \
#   -n no_predicted_70
# done

# for i in {1..4}; do
#   SAMPLE=$(sed -n "${i}p" dir_names.txt)
#   rm -r ${SAMPLE}/no_predicted_65
#   Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
#   -i ${SAMPLE}/annotation_no_predicted_locus_table_reps.txt \
#   -o ${SAMPLE} \
#   --allelic-bias 0.65 \
#   -n no_predicted_65
# done

# for i in {1..4}; do
#   SAMPLE=$(sed -n "${i}p" dir_names.txt)
#   rm -r ${SAMPLE}/predicted_70
#   Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
#   -i ${SAMPLE}/annotation_locus_table_reps.txt \
#   -o ${SAMPLE} \
#   --allelic-bias 0.7 \
#   -n predicted_70
# done

# for i in {1..4}; do
#   SAMPLE=$(sed -n "${i}p" dir_names.txt)
#   rm -r ${SAMPLE}/predicted_65
#   Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
#   -i ${SAMPLE}/annotation_locus_table_reps.txt \
#   -o ${SAMPLE} \
#   --allelic-bias 0.65 \
#   -n predicted_65
# done



