#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J sc_Allelome.PRO2.snRNAseq
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1   
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --array=1-25

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 9w/sc_filelist.txt)
echo "Processing sample: "$SAMPLE

annotation=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/chr_annotation_mm39.bed
SNPfile=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.bed

/usr/bin/time -v bash /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Allelome.LINK/00_src/Allelome.PRO2.sh \
  -i "$SAMPLE" \
  -a $annotation \
  -s $SNPfile \
  -o /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_heart_snRNAseq/9w/sc_Allelome.PRO2 \
  -r 1 \
  -t 10

exit_status=$?

echo
echo "===== SLURM ACCOUNTING SUMMARY ====="
jobref="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

sacct -M serial -r serial_std \
  --jobs="$jobref" \
  --format=JobID,JobName%20,State,ExitCode,Elapsed,TotalCPU,ReqCPUS,ReqMem,MaxRSS,AveRSS,MaxVMSize

echo
echo "===== SLURM EFFICIENCY SUMMARY ====="
seff "$jobref" || true

exit $exit_status