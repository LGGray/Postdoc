#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_spase
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=1-2

# Per-gene escape from the inactive X, using spASE's beta-binomial model.
# Array task 1 = 9w, task 2 = 78w.
#
# WHAT THIS ANSWERS THAT THE TILE ANALYSIS DOES NOT. The chrX aggregate is
# built against whole-chromosome intervals, so it has already summed over gene
# identity: 12.7% escape is one number for the whole chromosome. This job keeps
# gene identity and fits a beta-binomial per gene, which is the null
# NEXT_ANALYSIS task 7 says the data needs - the autosomal ratio here is 15-35x
# overdispersed relative to binomial, and spASE fits that overdispersion per
# gene rather than assuming it.
#
# TWO STEPS:
#   1. A counting pass with --gene-bin-out, which emits gene x pixel allelic
#      counts on a 16um grid. Nothing in the pipeline emitted this before:
#      --locus-out pools over the whole section and the main table sums over
#      gene identity.
#   2. spatial/spase_scase.R, which fits scase() per gene, calibrates the
#      B6-ward mapping bias on the autosomal genes, and tests each chrX gene
#      against the false-escape floor.
#   3. spatial/spase_spatial.R, which fits spase() - a 2D thin-plate spline per
#      gene - and asks whether a gene's escape VARIES across the section. Its
#      autosomal control set is depth-matched and is what calibrates the test;
#      if that control is not uniform, the script says so and refuses to plot
#      fitted surfaces. Expect a null result: rho at 64um was 0.006/0.003 and
#      C(d) is flat, so escape does not vary spatially in the aggregate.
#      Note that a hit here would be REGIONAL variation in escape, not an XCI
#      patch - this genotype has no mosaic.
#
# Stage 3 is much slower than stage 2 (a spline per gene, over every pixel, for
# the chrX genes AND three depth-matched controls each). Run 'scase' alone first
# if you only want the per-gene escape numbers.
#
# ON RUNNING ITS OWN COUNTING PASS. The UMI calls are identical to the sweep's:
# same SNP bed, same --snp-offset, --min-mapq, --min-baseq and count unit, and
# the counter is deterministic given those. What must NOT differ is the SNP bed
# or the offsets, which is why they are set from the same values here rather
# than being parameterised. If you would rather have one pass produce
# everything, add --gene-bin-out and --gene-bed to spatial_ase_sweep.slurm and
# force a recount there; this job will then reuse its output.
#
# WHY 16um PIXELS - AND WHY THAT REASONING WAS WRONG. The original argument was
# that scase() fits a free overdispersion per gene, which needs MANY pixels
# more than it needs depth per pixel, so 16um (thousands of pixels) beats 64um
# (a few hundred).
#
# That is backwards, and spase_se_diagnostic.R measures why. At 16um, 90% of
# occupied (gene, pixel) cells hold exactly ONE molecule (mean 1.12). A
# beta-binomial estimates overdispersion from variation between clusters at a
# given cluster size, and BetaBinom(1, a, b) is Bernoulli(p) for EVERY value of
# phi - so with clusters of size 1 phi drops out of the likelihood entirely. It
# is not identified. That is why half the fitted phis sit exactly on 0 or 1 and
# why the reported SE does not move with phi.
#
# The consequence is NOT that the CIs are wrong: a parametric bootstrap puts
# the true SE at 1.03x the reported one, so they are fine. The consequence is
# that scase() at 16um is an expensive binomial, and this analysis must not be
# described as accounting for overdispersion. Identifying overdispersion needs
# replication WITHIN a pixel, i.e. LARGER pixels. Run 32 or 64um to give the
# model something to fit, and compare.
#
# The pixel size is still not the tile size of the tile analysis, and the
# 15-35x overdispersion measured there is between-tile variance at a much
# coarser unit - a real, separate phenomenon that this grid cannot see.

# ---------------------------------------------------------------------------
# ONE-OFF: INSTALLING spASE. Do this before the first run.
#
# spASE's DESCRIPTION declares NO Imports, so install_github pulls in nothing
# and a missing dependency surfaces as "lazy loading failed for package spASE",
# naming ONE package at a time, rather than as a missing-package error up
# front. So install them first.
#
# The list is spASE's NAMESPACE - the packages that must resolve for the build
# to finish. Do NOT use the wider set that a grep of the source produces: that
# also catches Seurat, caret, Rfast, fields, pals, quadprog and CompQuadForm,
# which are referenced only inside functions this pipeline never calls (Seurat
# in a 10x .h5 reader and in a class() string test) or are spacexr's
# dependencies rather than spASE's. They are also by far the slowest to build,
# and installing them as one batch is what makes the install fail.
#
#   conda activate RNAseq
#   R
#   options(repos = c(CRAN = "https://cloud.r-project.org"))   # else a Tk
#   #   mirror-chooser window opens and the terminal appears to hang
#   need <- c("Matrix","aod","doParallel","doSNOW","dplyr","foreach","ggplot2",
#             "latex2exp","viridis","RColorBrewer","ggExtra","ggrepel","mgcv",
#             "reshape2")
#   install.packages(need[!need %in% rownames(installed.packages())])
#   devtools::install_github("lulizou/spASE", upgrade="never")
#
# spacexr is a separate install (devtools::install_github("dmcable/spacexr"),
# or conda, which is what worked here). It is loaded before spASE because the
# README says spASE overrides several of its functions - that matters for
# RCTD/C-SIDE, which this pipeline does not use.
#
# spacexr's own dependencies fields/Rfast/spam/dotCall64 DO need a Fortran
# compiler. Where there is no gfortran, install them from binary
# (options(pkgType="binary")) or from conda-forge as r-fields, r-rfast.
# ---------------------------------------------------------------------------

module load slurm_setup

# --export=NONE means variables set on the sbatch command line never reach the
# job. Positional arguments do:
#
#   sbatch spatial_spase.slurm                  # 16um pixels, no_Xist mask
#   sbatch spatial_spase.slurm 16 force         # redo the counting pass
#   sbatch spatial_spase.slurm 8                # 8um pixels
#   sbatch spatial_spase.slurm 16 "" no_Xic500kb   # the Xic-masked SNP bed
#   sbatch spatial_spase.slurm 16 "" "" scase      # per-gene escape only
#   sbatch spatial_spase.slurm 16 "" "" spatial    # the 2D test only
#   sbatch spatial_spase.slurm 16 "" "" both 8     # spatial test at df=8
#   sbatch spatial_spase.slurm 16 "" "" sediag     # are the SEs right at all?
#   sbatch spatial_spase.slurm 16 "" "" both 5 64  # merge to 64um, no recount
#
# ARGUMENT 6, REBIN. Merges the 16um table into coarser pixels in memory rather
# than recounting: a pixel index is array %/% (um/2), so (a %/% 8) %/% 4 is
# a %/% 32 - exactly the grid a native --gene-bin-um 64 pass would write.
# Lossless, and minutes instead of a BAM pass. Must be a multiple of arg 1.
#
# It changes phi, the CI width and the 2D spatial fit. It does NOT change the
# escape numbers: sum(alt)/sum(ref+alt) per gene is invariant to how pixels are
# grouped, and this is asserted at load time. Nor does it rescue the
# beta-binomial - cells holding one molecule are 90% at 16um and still 70% at
# 64um. Use it to line the spatial stage up with the 64um tile analysis, which
# is the one thing it genuinely buys.
GENE_BIN_UM="${1:-${GENE_BIN_UM:-16}}"
FORCE="${2:-${FORCE:-}}"
SNP_LABEL="${3:-${SNP_LABEL:-no_Xist}}"
# scase | spatial | both | sediag - which stages to run after the counting pass.
# sediag reads the scase table rather than producing one, so run scase first.
STAGE="${4:-${STAGE:-both}}"
# Spline degrees of freedom for the spatial stage. spASE's default is 5 and its
# docs say to try several; each df writes its own output files.
DF="${5:-${DF:-5}}"
# Coarser pixels without recounting. 0 or empty = use the table as counted.
REBIN="${6:-${REBIN_UM:-0}}"

if [ -n "$REBIN" ] && [ "$REBIN" != "0" ]; then
  if [ $((REBIN % GENE_BIN_UM)) -ne 0 ] || [ "$REBIN" -lt "$GENE_BIN_UM" ]; then
    echo "REBIN (arg 6, $REBIN) must be a multiple of the pixel size (arg 1, $GENE_BIN_UM)." >&2
    echo "      Pixels can only be merged, never split." >&2
    exit 1
  fi
fi

case "$STAGE" in
  scase|spatial|both|sediag) ;;
  *) echo "Unknown stage '$STAGE' (arg 4): expected scase, spatial, both or sediag" >&2
     exit 1 ;;
esac

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
REPO=$BASE/Postdoc

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
  SAMPLE=9w
else
  SAMPLE=78w
fi

# Mirrors the SUF convention in the other spatial scripts: the default mask
# keeps every path as it was, any other mask writes alongside it so two masks
# can be compared rather than one silently replacing the other.
if [ "$SNP_LABEL" = "no_Xist" ]; then SUF=""; else SUF="_$SNP_LABEL"; fi
SNPS=$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_${SNP_LABEL}.bed

OUT_DIR=$SPATIAL/ase/$SAMPLE
GENE_BINS=$OUT_DIR/gene_bins_${GENE_BIN_UM}um${SUF}.tsv.gz
BAM=$SPATIAL/$SAMPLE/outs/possorted_genome_bam.bam
SPATIAL_DIR=$SPATIAL/$SAMPLE/outs/binned_outputs/square_002um/spatial

# Gene bodies, one interval per gene, every chromosome, gene name in column 4.
# Written by OCM_heart/core_escape_SNPs.R, which is the one place the coordinate
# conventions are handled. If it is missing, rebuild the beds:
#   cd $BASE/GRCm39 && OUT_DIR=. Rscript $REPO/OCM_heart/core_escape_SNPs.R
GENE_BED=$BASE/GRCm39/annotation_us_mm39_gene_level.bed

mkdir -p "$OUT_DIR"
if [ -n "$REBIN" ] && [ "$REBIN" != "0" ]; then EFF_UM=$REBIN; else EFF_UM=$GENE_BIN_UM; fi
echo "Sample $SAMPLE | counted at ${GENE_BIN_UM}um | analysed at ${EFF_UM}um | SNP mask $SNP_LABEL"

for f in "$BAM" "$SNPS" "$GENE_BED"; do
  if [ ! -s "$f" ]; then
    echo "Missing required input: $f" >&2
    [ "$f" = "$GENE_BED" ] && echo \
      "      Rebuild it: cd $BASE/GRCm39 && OUT_DIR=. Rscript $REPO/OCM_heart/core_escape_SNPs.R" >&2
    exit 1
  fi
done

# Column 4 must hold gene names. A bed whose column 4 is empty or "." would
# produce one interval named chrom:start-end per gene, which the R script would
# then happily analyse - so check here rather than discover it in the results.
name_col=$(awk -F'\t' 'NR==1{print $4}' "$GENE_BED")
if [ -z "$name_col" ] || [ "$name_col" = "." ]; then
  echo "$GENE_BED has no gene name in column 4 (got '$name_col')" >&2
  exit 1
fi
echo "Gene bed: $GENE_BED  ($(wc -l < "$GENE_BED") genes, first name '$name_col')"

# ---------------------------------------------------------------------------
# 1. COUNTING PASS
# ---------------------------------------------------------------------------
if [ -s "$GENE_BINS" ] && [ "$FORCE" != "1" ] && [ "$FORCE" != "force" ]; then
  echo "Reusing existing $GENE_BINS ($(( $(zcat "$GENE_BINS" | wc -l) - 1 )) gene x pixel rows; pass 'force' as arg 2 to redo)"
else
  source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
  conda activate sinto     # pysam lives here alongside sinto itself

  # --snp-offset -1 is not a guess: the offset probe scored -1 at 100.0%
  # against +0 at 49.0% on 200k chrX reads. Column 2 of this SNP bed is a
  # 1-based VCF position and pysam reports 0-based.
  #
  # No --require-gene, deliberately. Molecules with no gene body stay in the
  # main table, which is chromosome-wide by design, and are simply absent from
  # the gene table. Requiring a gene here would change what the main table
  # measures for no benefit.
  #
  # --out goes to a separate filename from the sweep's bin_allele_counts.tsv:
  # the two are byte-identical in content but this job must never overwrite
  # the sweep's table or its provenance sidecar.
  python3 "$REPO/spatial/ase_bin_allele_counts.py" \
    --bam "$BAM" \
    --spatial-dir "$SPATIAL_DIR" \
    --snps "$SNPS" \
    --out "$OUT_DIR/bin_allele_counts_spase${SUF}.tsv" \
    --provenance-out "${GENE_BINS%.tsv.gz}.provenance.tsv" \
    --x-chrom chrX \
    --autosome-chroms chr19,chr18,chr17 \
    --snp-offset -1 \
    --min-mapq 255 \
    --min-baseq 20 \
    --gene-bin-out "$GENE_BINS" \
    --gene-bin-um "$GENE_BIN_UM" \
    --gene-bed "$GENE_BED" \
    --threads "$SLURM_CPUS_PER_TASK"
  rc=$?
  conda deactivate
  if [ $rc -ne 0 ]; then echo "Counting pass failed" >&2; exit $rc; fi
fi

# ---------------------------------------------------------------------------
# 2/3. THE FITS
# ---------------------------------------------------------------------------
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

# Both scripts source spatial/spase_common.R for config, input loading and the
# matrix build, so they cannot disagree about which allele matrix1 holds.
run_stage () {
  local script="$1"; shift
  echo
  echo "=== $script ==="
  env SAMPLE="$SAMPLE" BASE="$SPATIAL" ASE_DIR="$OUT_DIR" OUT_DIR="$OUT_DIR" \
    GENE_BINS="$GENE_BINS" GENE_BIN_UM="$GENE_BIN_UM" SNP_LABEL="$SNP_LABEL" \
    CORES="$SLURM_CPUS_PER_TASK" DF="$DF" REBIN_UM="$REBIN" "$@" \
    Rscript "$REPO/spatial/$script"
}

rc=0
if [ "$STAGE" = "scase" ] || [ "$STAGE" = "both" ]; then
  run_stage spase_scase.R || rc=$?
  if [ $rc -ne 0 ]; then
    conda deactivate; echo "spase_scase.R failed" >&2; exit $rc
  fi
fi

if [ "$STAGE" = "spatial" ] || [ "$STAGE" = "both" ]; then
  run_stage spase_spatial.R || rc=$?
  if [ $rc -ne 0 ]; then
    conda deactivate; echo "spase_spatial.R failed" >&2; exit $rc
  fi
fi

# Does the beta-binomial's reported standard error actually respond to the
# overdispersion it fitted? In the first run it did not - SE / binomial SE was
# 1.00-1.05 at every value of phi - which would make every CI and q-value in
# the scase table binomial. Not part of `both`: it reads the scase table rather
# than producing one, and it is a question about the estimator, asked once.
if [ "$STAGE" = "sediag" ]; then
  if [ ! -s "$OUT_DIR/spase_scase_${SAMPLE}_${EFF_UM}um${SUF}.tsv" ]; then
    conda deactivate
    echo "No scase table to diagnose. Run the scase stage first:" >&2
    echo "  sbatch $0 $GENE_BIN_UM \"\" $SNP_LABEL scase" >&2
    exit 1
  fi
  run_stage spase_se_diagnostic.R || rc=$?
  if [ $rc -ne 0 ]; then
    conda deactivate; echo "spase_se_diagnostic.R failed" >&2; exit $rc
  fi
fi
conda deactivate

echo
echo "Outputs in $OUT_DIR:"
ls -la "$OUT_DIR"/spase_scase_${SAMPLE}_${EFF_UM}um${SUF}.* \
       "$OUT_DIR"/spase_spatial_${SAMPLE}_${EFF_UM}um_df${DF}${SUF}.* \
       "$OUT_DIR"/spase_se_diagnostic_${SAMPLE}_${EFF_UM}um${SUF}.* \
       2>/dev/null
echo
echo "n = 1 animal per age: do not report a difference between the two"
echo "sections as an age effect."
