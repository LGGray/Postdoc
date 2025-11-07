#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J adult_Allelome_LINK
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
#SBATCH --array=1-27

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_bodymap

sample_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/01_aging_TH_for_LG/00_rawdata/01_adult
genome_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/STAR2.7.11b_genome_index/genome_dir_gencodeM36_mm39.202410
annotation=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/20250512_RefSeq.gtf
annotation_f=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation_f.bed
annotation_r=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation_r.bed
SNPfile=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed

run_one() {
  SAMPLE="$1"
  filename=$(echo "$SAMPLE" | cut -d'_' -f2)

  mkdir -p "adult/${SAMPLE}"

  shopt -s nullglob
  R1=($(printf "%s\n" ${adult_dir}/${SAMPLE}/${filename}_S*_L*_R1_001.fastq.gz | sort))
  R2=($(printf "%s\n" ${adult_dir}/${SAMPLE}/${filename}_S*_L*_R2_001.fastq.gz | sort))
  shopt -u nullglob
  IFS=, R1j="${R1[*]}"; unset IFS
  IFS=, R2j="${R2[*]}"; unset IFS

  # One STAR per srun; --exclusive ensures each gets its own CPU set
  srun --exclusive -n 1 -c "${SLURM_CPUS_PER_TASK}" \
    STAR \
      --genomeDir $genome_dir \
      --readFilesIn "$R1j" "$R2j" \
      --readFilesCommand gunzip -c \
      --runThreadN 12 \
      --outStd Log \
      --outSAMtype BAM SortedByCoordinate \
      --alignIntronMax 100000 \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outFileNamePrefix "adult/${SAMPLE}/${filename}_" \
      --outFilterMultimapNmax 1

  samtools index "adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam"

  STRAND=reverse
  htseq-count -s $STRAND -r pos -f bam \
    "adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam" \
    $annotation \
    > "adult/${SAMPLE}/${filename}_stranded_${STRAND}.count"

  echo "Finished processing sample: $i $SAMPLE"
}


SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" adult_path.txt)
filename=$(echo "$SAMPLE" | cut -d'_' -f2)
run_one "$SAMPLE"


# # Separate BAM files into forward and reverse reads
# BAM="adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam"
# ANNOTATION=../GRCm39/gencode.vM37.primary_assembly.annotation.gtf
# #awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' $ANNOTATION | gtf2bed - > ../GRCm39/gencode.vM37.primary_assembly.annotation.bed
# BED=../GRCm39/gencode.vM37.primary_assembly.annotation.bed

# ../Postdoc/infer_experiment.py -i $BAM -r $BED

strandrule="1+-,1-+,2++,2--"
perl /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/separate_BAM_strand.pl "adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out.bam" $strandrule adult/${SAMPLE}

# Split RefSeq into forward and reverse strands
#awk '$6 == "+" {print > "../GRCm39/forward_20250512_RefSeq.bed"} $6 == "-" {print > "../GRCm39/reverse_20250512_RefSeq.bed"}' ../GRCm39/20250512_RefSeq.bed

# Run Allelome.PRO
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
  -a $annotation_f \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o "adult/${SAMPLE}"

bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_rev.bam \
  -a $annotation_r \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o "adult/${SAMPLE}"

# Merge locus_table from forward and reverse strands
DATE=$(date +%Y_%m_%d)
Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R adult/$SAMPLE $DATE annotation


# # Allelome Pro with no predicted genes
# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
#   -a $ref_dir/annotation_no_predicted_f.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o "adult/${SAMPLE}"

# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_rev.bam \
#   -a $ref_dir/annotation_no_predicted_r.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o "adult/${SAMPLE}"

# # Allelome Pro with predicted genes
# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
#   -a $ref_dir/annotation_f.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o "adult/${SAMPLE}"

# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i adult/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_rev.bam \
#   -a $ref_dir/annotation_r.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o "adult/${SAMPLE}"

# # Merge locus_table from forward and reverse strands
# DATE=$(date +%Y_%m_%d)
# Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R adult/$SAMPLE $DATE annotation_no_predicted
# Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R adult/$SAMPLE $DATE annotation
