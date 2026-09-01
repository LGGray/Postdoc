#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_gene_escape_qc
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00

# QC layer for the per-gene chrX escape screen: SNP support, counting level,
# tile recurrence, between-age concordance. Complements spatial/escape_genes.R
# rather than replacing it - see the header of spatial/chrX_escape_qc.R.
#
#   sbatch ~/Postdoc/slurm/spatial_gene_escape_qc.slurm
#   sbatch ~/Postdoc/slurm/spatial_gene_escape_qc.slurm 64 pysam_dup
#   sbatch ~/Postdoc/slurm/spatial_gene_escape_qc.slurm 64 pysam        # dedup run
#
# Read-only: it reads tile_gene_counts.tsv[.gz] and the SNP/annotation beds, and
# writes only the three sidecars next to the input. Safe against a live run.
#
# 32G because tile_gene_counts.tsv is ~200M rows uncompressed for the dedup run
# and fread holds the chrX subset plus the 626k-row chrX SNP table. The dup run
# is gzipped and smaller. 30 min is generous; the foverlaps is the slow part and
# it is seconds.

module load slurm_setup

# Positional, not environment: --export=NONE stops `TILE_UM=64 sbatch` from
# reaching the job. Must match the tile run being read.
TILE_UM="${1:-64}"
# Which tile run to read. Names the input directory ase_<label>_<um>um, so it is
# the SNP_LABEL the pysam tile run used - "pysam_dup" for the duplicate-inclusive
# run, "pysam" for the UMI-collapsed one. NOTE the script reads the *_umi columns
# either way, so the two runs should agree closely; a large difference between
# them is a finding, not a formatting choice.
SNP_LABEL="${2:-pysam_dup}"
SAMPLES="${3:-9w,78w}"

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
REPO=$BASE/Postdoc
IN_DIR=$BASE/adult_aged_spatial/ase_${SNP_LABEL}_${TILE_UM}um

if [ ! -d "$IN_DIR" ]; then
  echo "ERROR: no such input directory: $IN_DIR" >&2
  echo "       available:" >&2
  ls -d "$BASE"/adult_aged_spatial/ase_* 2>/dev/null | sed 's/^/         /' >&2
  exit 1
fi

# The SNP bed is used only to COUNT informative SNPs per gene. It must be the bed
# the tile run was scored against, or the counts describe a different experiment.
# The pysam runs used the no_Xist bed; pass SNP_BED to override.
SNP_BED="${SNP_BED:-$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed}"

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, NOT the repo-wide RNAseq env. Most wrappers in this directory use
# RNAseq and are misleading here; same situation as spatial_tile_map.slurm.
# Do not "fix" this for consistency with the neighbours.
conda activate seurat_env

TILE_UM="$TILE_UM" SNP_LABEL="$SNP_LABEL" SAMPLES="$SAMPLES" \
  IN_DIR="$IN_DIR" SNP_BED="$SNP_BED" \
  ANNOT="$BASE/GRCm39/annotation_us_mm39_chrX.bed" \
  CORE_BED="$BASE/GRCm39/annotation_us_mm39_core_escape.bed" \
  Rscript "$REPO/spatial/chrX_escape_qc.R"

echo "Done. PDF, TSV and provenance in $IN_DIR"
