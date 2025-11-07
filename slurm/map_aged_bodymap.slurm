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
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --array=1-24


module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_bodymap

genome_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/STAR2.7.11b_genome_index/genome_dir_gencodeM36_mm39.202410
annotation=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/20250512_RefSeq.gtf
annotation_f=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation_f.bed
annotation_r=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation_r.bed
SNPfile=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed


run_one() {
  SAMPLE="$1"
  filename=$(basename "$SAMPLE" | cut -d'_' -f2)

  tmp=$(basename ${SAMPLE})
  mkdir -p "aged/${tmp}"

  # Build R1/R2 inputs: special-case samples need sorted, comma-joined lists
  shopt -s nullglob
  if [[ "$filename" == "22L012008" || "$filename" == "22L012009" || "$filename" == "22L012010" ]]; then
    R1=($(printf "%s\n" ${SAMPLE}/${filename}_S*_L*_R1_001.fastq.gz | sort))
    R2=($(printf "%s\n" ${SAMPLE}/${filename}_S*_L*_R2_001.fastq.gz | sort))
    IFS=, R1j="${R1[*]}"; unset IFS
    IFS=, R2j="${R2[*]}"; unset IFS

    STAR \
      --genomeDir $genome_dir \
      --readFilesIn "$R1j" "$R2j" \
      --readFilesCommand gunzip -c \
      --runThreadN 16 \
      --outStd Log \
      --outSAMtype BAM SortedByCoordinate \
      --alignIntronMax 100000 \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outFileNamePrefix "aged/${tmp}/${filename}_" \
      --outFilterMultimapNmax 1
  else
    R1j="${SAMPLE}/${filename}_S*_L*_R1_001.fastq.gz"
    R2j="${SAMPLE}/${filename}_S*_L*_R2_001.fastq.gz"

        STAR \
      --genomeDir $genome_dir \
      --readFilesIn $R1j $R2j \
      --readFilesCommand gunzip -c \
      --runThreadN 16 \
      --outStd Log \
      --outSAMtype BAM SortedByCoordinate \
      --alignIntronMax 100000 \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outFileNamePrefix "aged/${tmp}/${filename}_" \
      --outFilterMultimapNmax 1
  fi

  samtools index "aged/${tmp}/${filename}_Aligned.sortedByCoord.out.bam"

  STRAND=reverse
  htseq-count -s $STRAND -r pos -f bam \
    "aged/${tmp}/${filename}_Aligned.sortedByCoord.out.bam" \
    $annotation \
    > "aged/${tmp}/${filename}_stranded_${STRAND}.count"

  echo "Finished processing sample: $tmp"
}

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" aged_path.txt)
filename=$(basename "$SAMPLE" | cut -d'_' -f2)
run_one "$SAMPLE"

# # Separate BAM files into forward and reverse reads
strandrule="1+-,1-+,2++,2--"
perl /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/separate_BAM_strand.pl "aged/Sample_${filename}/${filename}_Aligned.sortedByCoord.out.bam" $strandrule aged/Sample_${filename}

# Run Allelome.PRO
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i aged/Sample_${filename}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
  -a $annotation_f \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o "aged/Sample_${filename}"

bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i aged/Sample_${filename}/${filename}_Aligned.sortedByCoord.out_rev.bam \
  -a $annotation_r \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o "aged/Sample_${filename}"

# Merge locus_table from forward and reverse strands
DATE=$(date +%Y_%m_%d)
Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R aged/Sample_${filename} $DATE annotation

# # Allelome Pro with no predicted genes
# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
#   -a $ref_dir/annotation_no_predicted_f.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o "aged/${SAMPLE}"

# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_rev.bam \
#   -a $ref_dir/annotation_no_predicted_r.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o "aged/${SAMPLE}"

# # Allelome Pro with predicted genes
# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_fwd.bam \
#   -a $ref_dir/annotation_f.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o "aged/${SAMPLE}"

# bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
#   -i aged/${SAMPLE}/${filename}_Aligned.sortedByCoord.out_rev.bam \
#   -a $ref_dir/annotation_r.bed \
#   -s /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed \
#   -o "aged/${SAMPLE}"


# # Merge locus_table from forward and reverse strands
# DATE=$(date +%Y_%m_%d)
# Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R aged/$SAMPLE $DATE annotation_no_predicted
# Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/merge_strands.R aged/$SAMPLE $DATE annotation