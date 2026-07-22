#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J index_igv_bams
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --export=NONE

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM
TABLE="$BASE/Allelic_ratio_results/top_ventricular_cardiomyocyte_per_sample.txt"
OUTDIR="$BASE/Allelic_ratio_results/igv_bams"
mkdir -p "$OUTDIR"

MANIFEST="$OUTDIR/manifest.txt"
> "$MANIFEST"

tail -n +2 "$TABLE" | while IFS=$'\t' read -r sample selection cell_barcode barcode_stripped total_reads allelic_ratio; do
  bam="$BASE/sinto/$sample/${cell_barcode}.bam"
  if [ ! -f "$bam" ]; then
    echo "WARNING: BAM not found: $bam" >&2
    continue
  fi

  samtools index "$bam"

  out_name="${sample}_${selection}_${barcode_stripped}"
  cp "$bam" "$OUTDIR/${out_name}.bam"
  cp "${bam}.bai" "$OUTDIR/${out_name}.bam.bai"
  echo -e "${sample}\t${selection}\t${barcode_stripped}\t${out_name}.bam" >> "$MANIFEST"
done

echo "Done. Indexed/copied BAMs are in $OUTDIR"
