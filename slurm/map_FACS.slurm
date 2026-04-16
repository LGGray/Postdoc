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
#SBATCH --array=2-32

data_path=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/00_raw_data/20230619_Project_699_Eckstein_Engelhardt_Stefan_XBxC_MCL_78w_XX_XY
# path=20240507_Project_1050_Eckstein_Engelhardt_Stefan_XBxC_MCL_TAC

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd aged_cardiac_RNAseq/

genome_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/STAR2.7.11b_genome_index/genome_dir_gencodeM36_mm39.202410
annotation=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation.gtf
annotation_us=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation_us.bed

SNPfile=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/mouse_SNP/SNPfile_C57BL6_NJ_CAST_EiJ.bed

run_one() {
  SAMPLE="$1"
  # One STAR per srun; --exclusive ensures each gets its own CPU set
  STAR \
    --genomeDir $genome_dir \
    --readFilesIn $data_path/Sample_$SAMPLE/${SAMPLE}_S*_L001_R1_001.fastq.gz,$data_path/Sample_$SAMPLE/${SAMPLE}_S*_L002_R1_001.fastq.gz $data_path/Sample_$SAMPLE/${SAMPLE}_S*_L001_R2_001.fastq.gz,$data_path/Sample_$SAMPLE/${SAMPLE}_S*_L002_R2_001.fastq.gz \
    --readFilesCommand zcat \
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

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample.file.txt )
DATE=$(date +%Y_%m_%d)
mkdir -p $SAMPLE
run_one "$SAMPLE"


# Run Allelome.PRO
bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i ${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam \
  -a $annotation_us \
  -s $SNPfile \
  -r 1 \
  -t 1 \
  -o ${SAMPLE}


Rscript /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.LINK.R \
    -i ${SAMPLE}/${DATE}_${SAMPLE}_Aligned.sortedByCoord.out_annotation_us.bed_1/locus_table.txt\
    -o ${SAMPLE} \
    --allelic-bias 0.65 \
    -r 20 \
    -w 100 \
    -n AL_65 \
    -d FALSE


/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/00_raw_data/20230619_Project_699_Eckstein_Engelhardt_Stefan_XBxC_MCL_78w_XX_XY/Sample_23L003319/23L003319_S1_L001_R1_001.fastq.gz