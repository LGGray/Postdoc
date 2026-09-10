#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J derive_marker_panels
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
#SBATCH --time=3:00:00

# Derive spatial marker panels from the snRNA-seq annotation, doublet-free.
#
#   sbatch slurm/derive_marker_panels.slurm            # 50 genes per set
#   sbatch slurm/derive_marker_panels.slurm out 30     # 30 genes, into OCM/out
#
# --export=NONE, so the arguments are positional. Memory is readRDS on the
# 5 GB object plus two sparse passes over the counts.

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate seurat_env

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
cd $BASE/OCM

OUT_DIR=${1:-marker_panels}
N_PER_SET=${2:-50}

export DOUBLET_FILE=Allelic_ratio_results/scDblFinder_per_cell.txt
if [ ! -f "$DOUBLET_FILE" ]; then
  echo "Missing $DOUBLET_FILE - run slurm/scDblFinder_rates.slurm first." >&2
  exit 1
fi

echo "OUT_DIR:   $PWD/$OUT_DIR"
echo "N_PER_SET: $N_PER_SET"

Rscript "$BASE/Postdoc/OCM_heart/derive_marker_panels.R" "$OUT_DIR" "$N_PER_SET"

echo
echo "Then annotate with the derived panels:"
echo "  sbatch --export=ALL,BIN=square_008um,MARKER_CSV=$BASE/OCM/$OUT_DIR/marker_panels_snRNA.csv \\"
echo "    slurm/spatial_cell_annotation.slurm"
