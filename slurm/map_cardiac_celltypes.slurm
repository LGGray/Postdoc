#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J cardiac_celltypes_Allelome_LINK
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
#SBATCH --array=1-4

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/cardiac_RNAseq

genome_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/STAR2.7.11b_genome_index/genome_dir_gencodeM36_mm39.202410
annotation=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation.gtf
annotation_us=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation_us.bed

SNPfile=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed

run_one() {
  SAMPLE="$1"
  # One STAR per srun; --exclusive ensures each gets its own CPU set
  STAR \
    --genomeDir $genome_dir \
    --readFilesIn fastq/${SAMPLE}_1.fastq fastq/${SAMPLE}_2.fastq \
    --runThreadN 12 \
    --outStd Log \
    --outSAMtype BAM SortedByCoordinate \
    --alignIntronMax 100000 \
    --outFilterIntronMotifs RemoveNoncanonical \
    --outFileNamePrefix "${SAMPLE}/${SAMPLE}_" \
    --outFilterMultimapNmax 1

  samtools index "${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

  STRAND=no
  htseq-count -s $STRAND -r pos -f bam \
    "${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam" \
    $annotation \
    > "${SAMPLE}/${SAMPLE}_stranded_${STRAND}.count"

  echo "Finished processing sample: $SAMPLE"
}

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_files.txt )
mkdir -p $SAMPLE
# run_one "$SAMPLE"

# # # Separate BAM files into forward and reverse reads
# BAM="${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
# # ANNOTATION=../GRCm39/gencode.vM37.primary_assembly.annotation.gtf
# # #awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' $ANNOTATION | gtf2bed - > ../GRCm39/gencode.vM37.primary_assembly.annotation.bed
# BED=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/20250512_RefSeq.bed

# ../Postdoc/infer_experiment.py -i $BAM -r $BED

# # strandrule="1+-,1-+,2++,2--"
# perl /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/separate_BAM_strand.pl "${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam" $strandrule ${SAMPLE}

# # Split RefSeq into forward and reverse strands
# #awk '$6 == "+" {print > "../GRCm39/forward_20250512_RefSeq.bed"} $6 == "-" {print > "../GRCm39/reverse_20250512_RefSeq.bed"}' ../GRCm39/20250512_RefSeq.bed

# Run Allelome.PRO
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK_mod/00_src/Allelome.PRO2.sh \
  -i ${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam \
  -a $annotation_us \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o ${SAMPLE}


for i in {1..8}; do
  dir=$(sed -n "${i}p" sample_dirs.txt)
  Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
    -i ${dir}/annotation_us_locus_table_reps.txt\
    -o ${dir} \
    --allelic-bias 0.65 \
    -r 20 \
    -w 100 \
    -n AL_65 \
    -d FALSE

done