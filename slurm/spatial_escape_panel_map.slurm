#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_escape_panel_map
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

# The core escape panel POOLED per 64 um tile, 9w beside 78w, allelic ratio
# above and depth below - spatial/tile_escape_panel_map.R.
#
# NO COUNTING HAPPENS HERE. The pysam run was scored against the genome-wide
# annotation_us.bed, so the tile x gene table already carries all eleven of
# new_core_escape_genes, and the no_Xist SNP mask drops only the Xist gene span,
# so Ftx and Jpx survive it. This is a read of an existing table: minutes, not
# the 4 h the counter took.
#
#   sbatch ~/Postdoc/slurm/spatial_escape_panel_map.slurm                  # dup tree, both variants
#   sbatch ~/Postdoc/slurm/spatial_escape_panel_map.slurm ase_pysam_64um   # umi tree
#   sbatch ~/Postdoc/slurm/spatial_escape_panel_map.slurm ase_pysam_32um 32
#   sbatch ~/Postdoc/slurm/spatial_escape_panel_map.slurm "" 64 full       # skip the noXic variant
#   sbatch ~/Postdoc/slurm/spatial_escape_panel_map.slurm "" 64 both 1    # no depth cutoff
#
# DEPTH CUTOFF (arg 4, default 4) - pooled MOLECULES before a tile is coloured,
# the same MIN_DEPTH knob tile_gene_ar_maps.R and tile_gene_ar_panel.R take.
# Be ready for what 4 does to this panel: it leaves 9 of 1069 tiles at 9w and 5
# of 797 at 78w, and no tile whose molecules were all CAST survives it, so the
# maps come out nearly empty. That is the depth, not a broken figure - the whole
# panel pools to a median of one molecule per tile. Pass 1 for the ungated map,
# which is the one to read alongside it.
#
# WHY THIS IS A JOB AND NOT JUST AN Rscript. The figure itself is small, but
# tile_geometry() reads the 2 um tissue_positions for both sections - millions
# of bins per section - and fits an lm over them to recover the tile side in
# pixel units, which builds a model matrix of that height. That is the memory
# driver, not the counts, and it is why this asks for 48G like
# spatial_tile_gene_ar.slurm rather than the few GB the tables suggest. On a
# login node it is the step that gets the job killed.
#
# TWO VARIANTS BY DEFAULT, and the second is not cosmetic. Ftx, Jpx and
# 5530601H04Rik are Xic-linked: with B6 carrying the Xist deletion they are
# expressed FROM the inactive CAST X rather than escaping from it, and they sit
# at AR 0.17-0.65 against Pbdc1 at 0.87. Keeping them drags the pooled panel AR
# down, which makes a monoallelic threshold look easier to clear than it is.
#   full   all 11 genes, matching new_core_escape_genes as written
#   noXic  the 8 conventional escapees - the conservative panel
# Arg 3 selects: "both" (default), "full", or "noXic".
#
# MOLECULES, not reads, in both variants. On this --keep-duplicates tree the
# read columns are duplicate-inclusive AND the duplication is allele-asymmetric
# (B6 ~17x against CAST ~10x on chrX at 9w), so a read-level pooled ratio biases
# toward monoallelic B6 - the exact artefact this figure would be read as
# evidence of. LEVEL=reads exists in the script and warns; it is not run here.
#
# WHAT THE FIGURE CAN AND CANNOT SUPPORT, so it is in the job log as well as the
# script header: a 5%-level monoallelic (Xi-loss) call needs ~9 pooled molecules
# at 9w and ~12 at 78w given the panel AR, and no tile on either section reaches
# either number - the deepest holds five. Coarsening to 256/320/512 um super-
# tiles needs no recount and still leaves nothing at BH q < 0.05. The page is a
# description of where the panel's molecules landed and which allele they
# carried, not a map of X loss.
#
# Positional, not environment: --export=NONE means `TREE=... sbatch ...` sets the
# variable in the login shell and the job never sees it.

module load slurm_setup

# ${1:-...} substitutes on EMPTY as well as unset, so `sbatch ... "" 64 full`
# reaches arg 3 without a second guard for the empty first argument.
TREE="${1:-ase_pysam_dup_64um}"
TILE_UM="${2:-64}"
WHICH="${3:-both}"                            # both | full | noXic
MIN_DEPTH="${4:-4}"                           # pooled molecules to colour a tile

case "$MIN_DEPTH" in
  ''|*[!0-9]*) echo "MIN_DEPTH (arg 4) must be a non-negative integer, got '$MIN_DEPTH'" >&2; exit 1 ;;
esac

case "$WHICH" in
  both|full|noXic) ;;
  *) echo "Unknown variant '$WHICH' (arg 3): expected both, full or noXic" >&2; exit 1 ;;
esac

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
REPO=$BASE/Postdoc
IN_ROOT=$SPATIAL/$TREE

echo "tree=$TREE  tile=${TILE_UM}um  variants=$WHICH  MIN_DEPTH=$MIN_DEPTH molecules"

# -s is not enough on its own for the uncompressed form: check both spellings,
# because spatial_tile_locus_map.slurm writes .gz and an interrupted gunzip
# leaves the plain name behind.
for s in 9w 78w; do
  f=$IN_ROOT/$s/tile_gene_counts.tsv.gz
  if [ ! -s "$f" ] && [ ! -s "${f%.gz}" ]; then
    echo "Missing $f - run spatial_tile_locus_map.slurm for $s first." >&2; exit 1
  fi
done

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, not RNAseq: the spatial R scripts are built against it, and
# tile_ratio_map.R and spatial_tile_gene_ar.slurm pin the same one.
conda activate seurat_env

# arrow is fatal here - there is no .tsv/.csv beside these parquets. png costs
# only the H&E underlay, patchwork only the 2x2 layout (the script falls back to
# gridExtra, then to one panel per page).
Rscript -e 'for (p in c("arrow","png","patchwork")) if (!requireNamespace(p, quietly=TRUE))
              message("note: R package \"", p, "\" not installed")' 2>&1 | sed 's/^/  /'

rc=0

if [ "$WHICH" = "both" ] || [ "$WHICH" = "full" ]; then
  echo
  echo "=== all 11 genes (new_core_escape_genes as written) ==="
  IN_ROOT="$IN_ROOT" TILE_UM="$TILE_UM" SAMPLES="9w,78w" LEVEL=umi DROP_XIC=FALSE \
    MIN_DEPTH="$MIN_DEPTH" \
    Rscript "$REPO/spatial/tile_escape_panel_map.R" || rc=$?
fi

if [ "$WHICH" = "both" ] || [ "$WHICH" = "noXic" ]; then
  echo
  echo "=== 8 conventional escapees (Ftx, Jpx, 5530601H04Rik dropped) ==="
  IN_ROOT="$IN_ROOT" TILE_UM="$TILE_UM" SAMPLES="9w,78w" LEVEL=umi DROP_XIC=TRUE \
    MIN_DEPTH="$MIN_DEPTH" \
    Rscript "$REPO/spatial/tile_escape_panel_map.R" || rc=$?
fi

conda deactivate

# Reported rather than swallowed: with two variants in one job a silent failure
# of the second leaves the first one's PDF sitting there looking like success.
if [ $rc -ne 0 ]; then
  echo "tile_escape_panel_map.R failed (rc=$rc) - see above" >&2
  exit $rc
fi

echo
echo "Done. Outputs in $IN_ROOT:"
ls -la "$IN_ROOT"/tile_escape_panel_map_*
