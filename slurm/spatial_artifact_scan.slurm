#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_artifact_scan
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

# Where the chrX CAST signal comes from, and a mask for the part of it that
# cannot be escape. Array task 1 = 9w, task 2 = 78w.
#
# WHY. Splitting the SAME molecules by whether they sit inside an annotated
# gene body takes the 12.7% chrX figure apart:
#
#                        9w      78w
#   chrX, whole        0.1272  0.1264
#   chrX, gene bodies  0.0612  0.0613
#   chrX, outside      0.4750  0.4981    <- biallelic
#   autosomes, whole   0.4821  0.4841
#
# Non-genic chrX molecules are 16% of chrX and carry 60% of its CAST signal, at
# the autosomal value. In this genotype the B6 X is active in every cell, so a
# chrX locus reading 50/50 - let alone 100% CAST - is not reporting allelic
# expression. The autosomal genic/non-genic split barely moves (0.482 vs
# 0.485), so this is chrX-specific rather than a property of the annotation.
#
# A 100kb scan of the window table the counting pass ALREADY wrote localises it
# to a handful of loci: chrX:11.5-11.6Mb alone (no annotated gene) is 5 ref vs
# 6201 alt at 9w, 42% of all chrX CAST molecules on its own.
#
# TWO MODES, because they cost very different amounts:
#
#   windows   Runs only spatial/ase_artifact_scan.R against the existing
#             ase/<sample>/chrX_windows_64um.tsv.gz. No BAM pass. Minutes.
#             Enough to rank the regions and to see the genic/non-genic split
#             with its autosomal control. START HERE.
#   count     A fresh counting pass with --snp-out (the per-SNP allele ledger,
#             which nothing has emitted before), then the same R script with
#             pass 2 enabled. Hours - it is a full BAM pass. This is what
#             produces the mask bed, because a mask has to name SNPs.
#
#   sbatch spatial_artifact_scan.slurm                 # windows mode
#   sbatch spatial_artifact_scan.slurm count           # full, writes the mask
#   sbatch spatial_artifact_scan.slurm count force     # recount even if present
#   sbatch spatial_artifact_scan.slurm windows "" no_Xic500kb   # other mask
#
# WHAT THE AUTOSOMES ARE FOR. They are biallelic by construction in this cross,
# so an autosomal SNP or window sitting near 0 or 1 cannot be biology. Their
# rate at a given threshold IS the false positive rate of the identical call
# made on chrX, measured on the same reads by the same counter. That is why the
# counting pass below keeps chr19/18/17 in --snp-out-chroms and why the mask
# file deliberately retains its autosomal rows: dropping them would hide what
# the threshold costs.

module load slurm_setup

MODE="${1:-${MODE:-windows}}"
FORCE="${2:-${FORCE:-}}"
SNP_LABEL="${3:-${SNP_LABEL:-no_Xist}}"

case "$MODE" in
  windows|count) ;;
  *) echo "Unknown mode '$MODE' (arg 1): expected windows or count" >&2
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

# Same SUF convention as every other spatial script: the default mask keeps
# paths as they were, any other mask writes alongside it.
if [ "$SNP_LABEL" = "no_Xist" ]; then SUF=""; else SUF="_$SNP_LABEL"; fi
SNPS=$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_${SNP_LABEL}.bed

OUT_DIR=$SPATIAL/ase/$SAMPLE
BAM=$SPATIAL/$SAMPLE/outs/possorted_genome_bam.bam
SPATIAL_DIR=$SPATIAL/$SAMPLE/outs/binned_outputs/square_002um/spatial
GENE_BED=$BASE/GRCm39/annotation_us_mm39_gene_level.bed

SNP_LEDGER=$OUT_DIR/snp_ledger${SUF}.tsv.gz
# The window table the earlier sweep already wrote. Count mode writes its own
# alongside it rather than over it, for the same reason the spASE job does:
# this job must never clobber the sweep's outputs or their provenance.
WINDOWS=$OUT_DIR/chrX_windows_64um${SUF}.tsv.gz
WINDOWS_NEW=$OUT_DIR/artifact_windows_all_64um${SUF}.tsv.gz

mkdir -p "$OUT_DIR"
echo "Sample $SAMPLE | mode $MODE | SNP mask $SNP_LABEL"

for f in "$SNPS" "$GENE_BED"; do
  if [ ! -s "$f" ]; then
    echo "Missing required input: $f" >&2
    [ "$f" = "$GENE_BED" ] && echo \
      "      Rebuild it: cd $BASE/GRCm39 && OUT_DIR=. Rscript $REPO/OCM_heart/core_escape_SNPs.R" >&2
    exit 1
  fi
done

# ---------------------------------------------------------------------------
# COUNTING PASS - only in count mode
# ---------------------------------------------------------------------------
if [ "$MODE" = "count" ]; then
  if [ ! -s "$BAM" ]; then echo "Missing BAM: $BAM" >&2; exit 1; fi
  if [ -s "$SNP_LEDGER" ] && [ "$FORCE" != "1" ] && [ "$FORCE" != "force" ]; then
    echo "Reusing existing $SNP_LEDGER (pass 'force' as arg 2 to redo)"
  else
    source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
    conda activate sinto    # pysam lives here alongside sinto itself

    # Identical UMI calls to every other pass: same SNP bed, same
    # --snp-offset -1 (the offset probe scored -1 at 100.0% against +0 at
    # 49.0% on 200k chrX reads), same --min-mapq/--min-baseq, same count unit.
    # Nothing here is parameterised, so the ledger cannot describe a different
    # set of molecules from the tables it is meant to explain.
    #
    # --window-chroms all, not just chrX: the autosomal windows are the control
    # that says what a skewed window looks like when the truth is known.
    python3 "$REPO/spatial/ase_bin_allele_counts.py" \
      --bam "$BAM" \
      --spatial-dir "$SPATIAL_DIR" \
      --snps "$SNPS" \
      --out "$OUT_DIR/bin_allele_counts_artifact${SUF}.tsv" \
      --provenance-out "${SNP_LEDGER%.tsv.gz}.provenance.tsv" \
      --x-chrom chrX \
      --autosome-chroms chr19,chr18,chr17 \
      --snp-offset -1 \
      --min-mapq 255 \
      --min-baseq 20 \
      --gene-bed "$GENE_BED" \
      --snp-out "$SNP_LEDGER" \
      --snp-out-chroms all \
      --window-out "$WINDOWS_NEW" \
      --window-size 100000 \
      --window-chroms all \
      --threads "$SLURM_CPUS_PER_TASK"
    rc=$?
    conda deactivate
    if [ $rc -ne 0 ]; then echo "Counting pass failed" >&2; exit $rc; fi
  fi
  # Prefer the pass we just ran: it covers the autosomes, the old one is chrX
  # only, and the autosomal arm is the whole calibration.
  [ -s "$WINDOWS_NEW" ] && WINDOWS=$WINDOWS_NEW
fi

if [ ! -s "$WINDOWS" ]; then
  echo "No window table at $WINDOWS" >&2
  echo "      Run this job in 'count' mode to produce one." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# THE SCAN
# ---------------------------------------------------------------------------
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate RNAseq

env SAMPLE="$SAMPLE" BASE="$SPATIAL" ASE_DIR="$OUT_DIR" OUT_DIR="$OUT_DIR" \
  SNP_LABEL="$SNP_LABEL" WINDOWS="$WINDOWS" SNP_LEDGER="$SNP_LEDGER" \
  GENE_BED="$GENE_BED" WINDOW_SIZE=100000 \
  Rscript "$REPO/spatial/ase_artifact_scan.R"
rc=$?
conda deactivate
if [ $rc -ne 0 ]; then echo "ase_artifact_scan.R failed" >&2; exit $rc; fi

echo
echo "Outputs in $OUT_DIR:"
ls -la "$OUT_DIR"/artifact_windows_${SAMPLE}${SUF}.tsv \
       "$OUT_DIR"/artifact_snps_${SAMPLE}${SUF}.* \
       "$OUT_DIR"/artifact_scan_${SAMPLE}${SUF}.* 2>/dev/null
echo
echo "The mask bed is artifact_snps_<sample>.bed. It is NOT applied anywhere"
echo "automatically - build the next SNP file from it, give that file its own"
echo "SNP_LABEL, and re-run the sweep and spASE against both so the two can be"
echo "compared rather than one silently replacing the other."
