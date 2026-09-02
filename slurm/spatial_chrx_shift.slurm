#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_chrx_shift
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
#SBATCH --mem=24G
#SBATCH --time=02:00:00

# Which X-linked genes' allelic ratio differs between the two sections, and is
# each one's escape spatially patchy or scattered. Three steps:
#   spatial/chrx_allelic_shift.py       the scan, calibrated on the autosomes
#   spatial/chrx_escape_clustering.py   join-count permutation per hit gene
#   spatial/chrx_shift_figures.R        five pages, diagnostic first
#
# NO COUNTING HAPPENS HERE. This reads the tile x gene table that
# spatial_tile_locus_map.slurm already wrote, so it is minutes, not the 4 h the
# counter took. Memory is modest because nothing needs tile geometry - unlike
# spatial_tile_gene_ar.slurm, whose 48G goes on an lm over millions of 2 um bins.
#
#   sbatch ~/Postdoc/slurm/spatial_chrx_shift.slurm                    # dup tree, umi
#   sbatch ~/Postdoc/slurm/spatial_chrx_shift.slurm ase_pysam_64um     # umi tree
#   sbatch ~/Postdoc/slurm/spatial_chrx_shift.slurm "" umi 100         # looser depth floor
#   sbatch ~/Postdoc/slurm/spatial_chrx_shift.slurm "" umi 200 chr7    # imprinted loci instead
#
# READ THE CALIBRATION LINE BEFORE THE HIT LIST. Step 1 runs its own test on
# every autosomal gene as a null - autosomes cannot shift allelically between
# two sections, so whatever fraction of them comes back significant is this
# comparison's false-positive rate. Raw, it is ~7%. The script measures the
# overdispersion per depth stratum, divides it out, and prints the autosomal
# rate again; it should land near zero. If the corrected rate is still above a
# couple of percent the job says so, and nothing in the hit list is safe.
#
# LEVEL (arg 2, default umi). Molecules. --level reads is duplicate-inclusive
# and the duplication is allele-asymmetric, so it biases the ratio toward the
# more amplified allele - the script warns and runs it anyway, for comparison.
#
# DEPTH FLOOR (arg 3, default 200 molecules in BOTH sections). This is what
# excludes Sh3kbp1, Kctd12b and Plp1, which move on the escape panel but carry
# 138-157 molecules. Lower it and the inflation correction has to work harder;
# the calibration line is how you check whether it did.
#
# ONE SECTION PER AGE. n = 1 animal per group, so every difference here
# describes two sections and cannot be attributed to age. The scripts label
# their output with the sample names for that reason and never with
# adult/aged as groups.

module load slurm_setup

TREE="${1:-ase_pysam_dup_64um}"
LEVEL="${2:-umi}"
MIN_MOL="${3:-200}"
CHROM="${4:-chrX}"
EXTRA_ROOT="${5:-}"

case "$LEVEL" in
  umi|reads) ;;
  *) echo "LEVEL (arg 2) must be umi or reads, got '$LEVEL'" >&2; exit 1 ;;
esac
case "$MIN_MOL" in
  ''|*[!0-9]*) echo "MIN_MOL (arg 3) must be a non-negative integer, got '$MIN_MOL'" >&2; exit 1 ;;
esac

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
REPO=$BASE/Postdoc
IN_ROOT=$SPATIAL/$TREE

echo "tree=$TREE  level=$LEVEL  min_molecules=$MIN_MOL  chrom=$CHROM"

for s in 9w 78w; do
  f=$IN_ROOT/$s/tile_gene_counts.tsv.gz
  if [ ! -s "$f" ] && [ ! -s "${f%.gz}" ]; then
    echo "Missing $f - run spatial_tile_locus_map.slurm for $s first." >&2; exit 1
  fi
done

EXTRA_ARGS=()
if [ -n "$EXTRA_ROOT" ]; then
  for s in 9w 78w; do
    f=$EXTRA_ROOT/$s/tile_gene_counts.tsv.gz
    if [ ! -s "$f" ] && [ ! -s "${f%.gz}" ]; then
      echo "EXTRA_ROOT given but $f is missing." >&2; exit 1
    fi
  done
  EXTRA_ARGS=(--extra-root "$EXTRA_ROOT")
  echo "merging extra counts from $EXTRA_ROOT"
fi

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

# Step 1 and 2 are stdlib-only python - no pysam, no numpy - so the base
# interpreter is enough and the sinto env is not needed.
echo
echo "=== 1/3  scan, calibrated on the autosomal null ==="
python3 "$REPO/spatial/chrx_allelic_shift.py" \
  --in-root "$IN_ROOT" \
  --samples 9w,78w \
  --chrom "$CHROM" \
  --level "$LEVEL" \
  --min-molecules "$MIN_MOL" \
  "${EXTRA_ARGS[@]}"
rc=$?
if [ $rc -ne 0 ]; then echo "chrx_allelic_shift.py failed (rc=$rc)" >&2; exit $rc; fi

echo
echo "=== 2/3  spatial clustering of escape, per hit gene ==="
python3 "$REPO/spatial/chrx_escape_clustering.py" \
  --tiles "$IN_ROOT/chrx_shift_tiles_${LEVEL}.tsv.gz" \
  --shift-table "$IN_ROOT/chrx_allelic_shift_${LEVEL}.tsv" \
  --out-dir "$IN_ROOT"
rc=$?
if [ $rc -ne 0 ]; then echo "chrx_escape_clustering.py failed (rc=$rc)" >&2; exit $rc; fi

# seurat_env, not RNAseq: every spatial R script in this repo pins it.
conda activate seurat_env
Rscript -e 'if (!requireNamespace("ggrepel", quietly=TRUE))
              message("note: R package \"ggrepel\" not installed - labels may overlap")' 2>&1 | sed 's/^/  /'

echo
echo "=== 3/3  figures ==="
IN_ROOT="$IN_ROOT" LEVEL="$LEVEL" CHROM="$CHROM" \
  Rscript "$REPO/spatial/chrx_shift_figures.R"
rc=$?
conda deactivate
if [ $rc -ne 0 ]; then echo "chrx_shift_figures.R failed (rc=$rc)" >&2; exit $rc; fi

echo
echo "Done. Outputs in $IN_ROOT:"
ls -la "$IN_ROOT"/chrx_allelic_shift_* "$IN_ROOT"/chrx_escape_clustering* \
       "$IN_ROOT"/chrx_shift_tiles_* "$IN_ROOT"/chrx_shift_figures_* 2>/dev/null
