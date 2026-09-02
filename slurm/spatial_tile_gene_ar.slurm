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
#SBATCH --time=02:00:00

# Per-tile allelic ratio for a NAMED gene panel, 9w beside 78w. Two figures:
#   tile_gene_ar_panel.R   distributions - one violin/dot column per gene per age
#   tile_gene_ar_maps.R    slide maps    - each tile coloured by that gene's ratio
# Both read the same gene panel from spatial/tile_gene_panel.R, so they cannot
# end up showing different gene sets.
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
# DEPTH CUTOFF (arg 3, default 4) - MIN_DEPTH, units of the level being plotted.
# It stops a single molecule from painting a tile dark blue on the maps, and on
# the DISTRIBUTIONS it is a floor under MIN_TILE_N rather than a replacement, so
# it can only tighten: the read page stays at its existing 24, the umi page
# moves 2 -> 4.
#
# ON THE MAPS A READ CUTOFF OF 4 ONLY HALF WORKS, and it is worth knowing which
# half. A read here is a PCR duplicate and duplicates of one molecule all carry
# the SAME allele, so a tile with one molecule amplified four times passes. 53%
# of single-molecule tiles at 9w (60% at 78w) already have four or more reads,
# so MIN_DEPTH=4 at the read level keeps 52% of the CAST-only tiles. The cutoff
# that removes them is a MOLECULE cutoff - at the umi level, 4 keeps 0.7% of
# them - which is what arg 4 is for:
#   sbatch ~/Postdoc/slurm/spatial_tile_gene_ar.slurm "" 64 4  reads,umi
#   sbatch ~/Postdoc/slurm/spatial_tile_gene_ar.slurm "" 64 24 reads   # ~2 molecules
#
# MAP LEVELS (arg 4, default reads) - which levels get a map PDF. Each is its
# own file, named with the level and the cutoff, so nothing overwrites anything.
#
# Small enough to run interactively instead, which is the faster way to change
# the gene list:
#   conda activate seurat_env
#   Rscript ~/Postdoc/spatial/tile_gene_ar_panel.R
#
# Outputs, into the input tree so they sit beside the map they came from:
#   tile_gene_ar_panel.pdf              9 pages: {reads, umi} x {skew, imprinted,
#                                       escape}, then 3 tile-composition pages
#   tile_gene_ar_panel_pooled.csv       one row per gene per sample
#   tile_gene_ar_panel_per_tile.csv.gz  every tile behind the dots
#   tile_gene_ar_maps_64um_reads.pdf    one page per gene: allelic ratio across
#                                       the top, coverage underneath, ages left
#                                       to right, on the H&E
#
# The maps need arrow (tissue_positions.parquet) and png (the H&E underlay).
# Both are checked below; a missing png costs the H&E, a missing arrow is fatal
# because there is no .tsv/.csv fallback beside these parquets.

module load slurm_setup

TREE="${1:-ase_pysam_dup_64um}"
TILE_UM="${2:-64}"
MIN_DEPTH="${3:-4}"
MAP_LEVELS="${4:-reads}"

case "$MIN_DEPTH" in
  ''|*[!0-9]*) echo "MIN_DEPTH (arg 3) must be a non-negative integer, got '$MIN_DEPTH'" >&2; exit 1 ;;
esac

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

Rscript -e 'for (p in c("arrow","png","patchwork")) if (!requireNamespace(p, quietly=TRUE))
              message("note: R package \"", p, "\" not installed")' 2>&1 | sed 's/^/  /'

echo "=== distributions (MIN_DEPTH=$MIN_DEPTH, floored against MIN_TILE_N) ==="
IN_ROOT="$IN_ROOT" TILE_UM="$TILE_UM" SAMPLES="9w,78w" MIN_DEPTH="$MIN_DEPTH" \
  Rscript "$REPO/spatial/tile_gene_ar_panel.R"
rc=$?
if [ $rc -ne 0 ]; then conda deactivate; echo "tile_gene_ar_panel.R failed" >&2; exit $rc; fi

# One PDF per requested level. Before the cutoff existed this was one run at
# LEVEL=reads, on the grounds that a tile's colour is which allele its molecules
# came from and duplicating a molecule cannot change that. That reasoning is
# what BREAKS under a cutoff: the two levels now differ in which tiles are drawn
# at all, so whichever level the cutoff is meant to apply to has to be asked for.
IFS=',' read -r -a LEVS <<< "$MAP_LEVELS"
for lev in "${LEVS[@]}"; do
  lev="$(echo "$lev" | tr -d '[:space:]')"
  [ -z "$lev" ] && continue
  case "$lev" in
    reads|umi) ;;
    *) echo "Unknown map level '$lev' (arg 4): expected reads or umi" >&2; rc=1; continue ;;
  esac
  echo
  echo "=== slide maps: LEVEL=$lev MIN_DEPTH=$MIN_DEPTH ==="
  IN_ROOT="$IN_ROOT" TILE_UM="$TILE_UM" SAMPLES="9w,78w" LEVEL="$lev" \
    MIN_DEPTH="$MIN_DEPTH" Rscript "$REPO/spatial/tile_gene_ar_maps.R" || rc=$?
done

conda deactivate
if [ $rc -ne 0 ]; then echo "a maps or distributions run failed (rc=$rc)" >&2; exit $rc; fi

echo
echo "Done. Outputs in $IN_ROOT:"
ls -la "$IN_ROOT"/tile_gene_ar_*
