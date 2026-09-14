#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_tile_structure
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

# Quantify spatial structure in the per-tile chrX allelic ratio:
# global Moran's I with a permutation null, the same statistic after
# residualising on depth, LISA, and the beta-binomial tile null of
# NEXT_ANALYSIS task 7. See spatial/tile_spatial_structure.R for what each one
# can and cannot say - in particular that none of them is a test for XCI
# patches, because the CAST X is inactive in every cell of this cross.
#
# SIZING. serial_std caps cpus-per-task at 16 and sbatch REJECTS an
# over-request rather than trimming it (see the LRZ partition table). This job
# is single-threaded R over ~8000 tiles; 4 cpus and 32G is generous, and the
# 999-permutation run takes a couple of minutes, not hours. The 2h wall is
# slack for a larger NPERM or a finer tile size, not a measured need.
#
# Positional, not environment: --export=NONE stops `TREE=... sbatch` from
# reaching the job.
#
#   sbatch slurm/spatial_tile_structure.slurm
#   sbatch slurm/spatial_tile_structure.slurm ase_pysam_dup_64um read 20
#   sbatch slurm/spatial_tile_structure.slurm ase_pysam_64um umi 10 999
#
# Arg 1  TREE        counting tree under adult_aged_spatial/  (ase_pysam_dup_64um)
# Arg 2  COUNT_UNIT  read | umi                               (read)
# Arg 3  MIN_X       minimum chrX units per tile              (20)
# Arg 4  NPERM       permutations for the global null         (999)
# Arg 5  SAMPLES     comma separated                          (9w,78w)
#
# DUP_FACTOR is derived from COUNT_UNIT, not passed: a umi tree is already at
# the independent unit and needs no deflation, a read tree does. It only enters
# the beta-binomial null; Moran's I never references n.

module load slurm_setup

TREE="${1:-ase_pysam_dup_64um}"
COUNT_UNIT="${2:-read}"
MIN_X="${3:-20}"
NPERM="${4:-999}"
SAMPLES="${5:-9w,78w}"

case "$COUNT_UNIT" in
  read) DUP_FACTOR=3 ;;
  umi)  DUP_FACTOR=1 ;;
  *)    echo "COUNT_UNIT must be read or umi, got '$COUNT_UNIT'" >&2; exit 1 ;;
esac

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
REPO=$HOME/Postdoc
SPATIAL_BASE="$BASE/adult_aged_spatial"

if [ ! -d "$SPATIAL_BASE/$TREE" ]; then
  echo "no such tree: $SPATIAL_BASE/$TREE" >&2
  echo "available:" >&2; ls -d "$SPATIAL_BASE"/ase* >&2
  exit 1
fi

# The conda profile path is hardcoded, as in every other script here: under
# --export=NONE conda is not on PATH, so `conda info --base` cannot be used to
# find it.
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, NOT the repo-wide RNAseq env: this needs the same data.table and
# ggplot2 the rest of spatial/ uses, and the Matrix 1.8.0 dev tree on .libPaths.
conda activate seurat_env

echo "tree       $TREE"
echo "unit       $COUNT_UNIT  (dup factor $DUP_FACTOR)"
echo "min chrX   $MIN_X"
echo "perms      $NPERM"
echo "samples    $SAMPLES"

BASE="$SPATIAL_BASE" TREE="$TREE" SAMPLES="$SAMPLES" \
  COUNT_UNIT="$COUNT_UNIT" MIN_X="$MIN_X" NPERM="$NPERM" \
  DUP_FACTOR="$DUP_FACTOR" \
  TAG="${COUNT_UNIT}_min${MIN_X}" \
  Rscript "$REPO/spatial/tile_spatial_structure.R"

echo "Done. CSVs and PDF in $SPATIAL_BASE/$TREE/"
