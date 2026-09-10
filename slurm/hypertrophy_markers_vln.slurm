#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J hypertrophy_markers_vln
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
#SBATCH --time=2:00:00

# Hypertrophic marker violins for slide 17 of the IPT seminar deck: Sham vs
# TAC ventricular cardiomyocytes, deck palette, per-gene Wilcoxon rank-sum,
# scDblFinder doublets excluded.
#
#   sbatch slurm/hypertrophy_markers_vln.slurm              # -> OCM/hypertrophy_markers/
#   sbatch slurm/hypertrophy_markers_vln.slurm some_dir     # -> OCM/some_dir/
#
# Memory is dominated by readRDS on the 5 GB Seurat object; the plotting and
# the tests themselves are cheap. --export=NONE, so the argument is positional.

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, not the RNAseq env most of slurm/ uses: this is a Seurat script.
conda activate seurat_env

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
cd $BASE/OCM

OUT_DIR=${1:-hypertrophy_markers}

# The scDblFinder per-nucleus calls, written by OCM_heart/scDblFinder_rates.R.
export DOUBLET_FILE=Allelic_ratio_results/scDblFinder_per_cell.txt

if [ ! -f "$DOUBLET_FILE" ]; then
  echo "Missing $DOUBLET_FILE - run slurm/scDblFinder_rates.slurm first." >&2
  exit 1
fi
head -1 "$DOUBLET_FILE" | grep -q "class" || {
  echo "$DOUBLET_FILE has no 'class' column" >&2; exit 1; }

echo "DOUBLET_FILE: $DOUBLET_FILE"
echo "OUT_DIR:      $PWD/$OUT_DIR"

Rscript "$BASE/Postdoc/OCM_heart/hypertrophy_markers_vln.R" "$OUT_DIR"

echo "Done. Figure and stats in $PWD/$OUT_DIR"
