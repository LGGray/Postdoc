#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH -J adult_Allelome_LINK
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6          
#SBATCH --cpus-per-task=9
#SBATCH --mem=36G
#SBATCH --time=24:00:00

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_bodymap

adult_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/01_aging_TH_for_LG/00_rawdata/01_adult

# Read the first 27 samples from the list
mapfile -t SAMPLES < <(head -n 27 adult_path.txt)

run_one() {
  SAMPLE="$1"
  filename=$(echo "$SAMPLE" | cut -d'_' -f2)

  mkdir -p "adult/${SAMPLE}/HTseq_counts"

  shopt -s nullglob
  R1=($(printf "%s\n" ${adult_dir}/${SAMPLE}/${filename}_S*_L*_R1_001.fastq.gz | sort))
  R2=($(printf "%s\n" ${adult_dir}/${SAMPLE}/${filename}_S*_L*_R2_001.fastq.gz | sort))
  shopt -u nullglob
  IFS=, R1j="${R1[*]}"; unset IFS
  IFS=, R2j="${R2[*]}"; unset IFS

  # One STAR per srun; --exclusive ensures each gets its own CPU set
  srun --exclusive -n 1 -c "${SLURM_CPUS_PER_TASK}" \
    STAR \
      --genomeDir ../GRCm39/STAR_index \
      --readFilesIn "$R1j" "$R2j" \
      --readFilesCommand zcat \
      --runThreadN "${SLURM_CPUS_PER_TASK}" \
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
    ../GRCm39/gencode.vM37.primary_assembly.annotation.gtf \
    > "adult/${SAMPLE}/HTseq_counts/${filename}_stranded_${STRAND}.count"
}

# Launch up to --ntasks-per-node jobs concurrently
running=0
for S in "${SAMPLES[@]}"; do
  run_one "$S" &
  (( running++ ))
  if (( running >= SLURM_NTASKS_PER_NODE )); then
    wait
    running=0
  fi
done
wait