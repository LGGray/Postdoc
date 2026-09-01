#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_tile_gene_ar
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
#SBATCH --mem=48G
#SBATCH --time=01:00:00

# Per-tile allelic ratio for a NAMED gene panel, 9w beside 78w.
#
# Reads the tile x gene table spatial_tile_locus_map.slurm already wrote - no
# counting happens here, so this is minutes, not the 4 h the counter took. The
# default tree is the DUPLICATE-INCLUSIVE 64 um one:
#     adult_aged_spatial/ase_pysam_dup_64um/<sample>/tile_gene_counts.tsv.gz
# and that same file carries the deduplicated columns, so the umi page comes
# out of the same read rather than a second tree.
#
#   sbatch ~/Postdoc/slurm/spatial_tile_gene_ar.slurm                    # dup 64 um
#   sbatch ~/Postdoc/slurm/spatial_tile_gene_ar.slurm ase_pysam_64um     # umi tree
#   sbatch ~/Postdoc/slurm/spatial_tile_gene_ar.slurm ase_pysam_32um 32  # 32 um
#
# Small enough to run interactively instead, which is the faster way to change
# the gene list:
#   conda activate seurat_env
#   Rscript ~/Postdoc/spatial/tile_gene_ar_panel.R
#
# Outputs, into the input tree so they sit beside the map they came from:
#   tile_gene_ar_panel.pdf              6 pages: {reads, umi} x {skew, imprinted, escape}
#   tile_gene_ar_panel_pooled.csv       one row per gene per sample
#   tile_gene_ar_panel_per_tile.csv.gz  every tile behind the dots

module load slurm_setup

TREE="${1:-ase_pysam_dup_64um}"
TILE_UM="${2:-64}"

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
REPO=$BASE/Postdoc
IN_ROOT=$SPATIAL/$TREE

for s in 9w 78w; do
  f=$IN_ROOT/$s/tile_gene_counts.tsv.gz
  if [ ! -s "$f" ] && [ ! -s "${f%.gz}" ]; then
    echo "Missing $f - run spatial_tile_locus_map.slurm for $s first." >&2; exit 1
  fi
done

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, not RNAseq: the spatial R scripts are built against it, and
# tile_ratio_map.R already pins the same one.
conda activate seurat_env

IN_ROOT="$IN_ROOT" TILE_UM="$TILE_UM" SAMPLES="9w,78w" \
  Rscript "$REPO/spatial/tile_gene_ar_panel.R"
rc=$?
conda deactivate
if [ $rc -ne 0 ]; then echo "tile_gene_ar_panel.R failed" >&2; exit $rc; fi

echo
echo "Done. Outputs in $IN_ROOT:"
ls -la "$IN_ROOT"/tile_gene_ar_panel*
