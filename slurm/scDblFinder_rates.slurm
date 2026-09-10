#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J scDblFinder_rates
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=4:00:00

# Per-sample scDblFinder doublet rate on heart_seurat_object_SCT.rds.
# Read-only: reports the rate, does not rewrite the object or the sinto
# barcode lists. Re-run Seurat_preprocessing.R to actually drop the doublets.
#
# Memory is dominated by readRDS on the 5 GB object; scDblFinder itself works
# on one sample's sparse counts per worker.

module load slurm_setup

# --export=NONE means DBL_THREADS=... on the sbatch line does not reach the
# job. Pass the worker count as the first positional argument instead:
#
#   sbatch slurm/scDblFinder_rates.slurm        # 4 workers (default)
#   sbatch slurm/scDblFinder_rates.slurm 8      # 8 workers
#
# Keep it <= --cpus-per-task above, and resize both together.
DBL_THREADS="${1:-4}"

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, not the RNAseq env most of slurm/ uses: this is a Seurat script.
conda activate seurat_env

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
cd $BASE/OCM

if [ ! -f heart_seurat_object_SCT.rds ]; then
  echo "Missing heart_seurat_object_SCT.rds in $PWD" >&2
  exit 1
fi

DBL_THREADS="$DBL_THREADS" Rscript $BASE/Postdoc/OCM_heart/scDblFinder_rates.R
