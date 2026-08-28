#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_ase_sweep
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

# Step 0/1 of the spatial allele-specific analysis: measure how much allelic
# signal there actually is per unit area, and use it to pick a tile size.
#
# No BAM splitting happens here. Allelic counts are made once at the 2um grid,
# and every candidate tile size is derived from them by summation. sinto gets
# run afterwards, once, at the chosen size (spatial_sinto_tiles.slurm).
#
# Array task 1 = 9w, task 2 = 78w.
#
# Re-running is cheap: the counting pass is skipped if its output already
# exists AND is complete, so tuning the sweep is just `sbatch` again; pass
# 'force' to redo it. "Complete" now includes the window table and every
# requested subset column, so the first run after those were added recounts.
#
# Outputs, in adult_aged_spatial/ase/<sample>/:
#   bin_allele_counts.tsv                per-2um-bin counts, plus one column
#                                        pair per --subset-bed
#   bin_allele_counts.tsv.provenance.tsv SNP bed path and md5, informative SNP
#                                        and gene counts - written every run
#   chrX_windows_64um.tsv.gz             allele-split counts per 100 kb window
#                                        per 64um tile, for task 1
#   subset_locus_counts.tsv              per-interval UMI counts and ratios for
#                                        every named subset interval; this is
#                                        where you check the imprinted control
#                                        has depth, and which parent is B6
#   tile_size_sweep.{csv,pdf}            the sweep itself, now with a C(d) curve
#                                        per subset
#   pair_correlation_short_range.csv     C(4um) per set with the implied per-UMI
#                                        error rate - tasks 3 and 4 in one table
#   provenance_tile_sweep.tsv            what produced all of the above
#
# Then, once both array tasks are done:
#   sbatch ~/Postdoc/slurm/spatial_window_profile.slurm 64

module load slurm_setup

# --export=NONE above means variables set on the sbatch command line never
# reach the job - `TILE_UM=64 sbatch ...` sets it in your login shell and the
# job runs without it. Positional arguments do get through:
#
#   sbatch spatial_ase_sweep.slurm            # sweep only
#   sbatch spatial_ase_sweep.slurm 64         # also write the 64um sinto map
#   sbatch spatial_ase_sweep.slurm 64 force   # ...and redo the counting pass
#   sbatch spatial_ase_sweep.slurm 64 force 32   # ...window table on 32um tiles
#
# The environment forms still work if you submit with --export=ALL.
TILE_UM="${1:-${TILE_UM:-}}"
FORCE="${2:-${FORCE:-}}"
# Tile grid the window table is written on. Must match the TILE_UM that
# ase_chrX_window_profile.R is run at, or the join there finds nothing - it
# checks and says so rather than plotting an empty panel.
WINDOW_TILE_UM="${3:-64}"

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
REPO=$BASE/Postdoc

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
  SAMPLE=9w
elif [ "$SLURM_ARRAY_TASK_ID" -eq 2 ]; then
  SAMPLE=78w
else
  echo "Unexpected SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2; exit 1
fi

OUT_DIR=$SPATIAL/ase/$SAMPLE
COUNTS=$OUT_DIR/bin_allele_counts.tsv
mkdir -p "$OUT_DIR"

# Standard 10x GRCm39 reference, matching the rest of the lab's allelic work.
# This carries a small B6-ward mapping bias; the autosomal control set in the
# counter is there to measure it rather than to pretend it is absent.
BAM=$SPATIAL/$SAMPLE/outs/possorted_genome_bam.bam
SPATIAL_DIR=$SPATIAL/$SAMPLE/outs/binned_outputs/square_002um/spatial
SNPS=$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed

# ---- interval beds for the subset columns --------------------------------
#
# All written by OCM_heart/core_escape_SNPs.R, which is the one place the
# coordinate conventions are handled. Run it once before this job:
#   cd $BASE/GRCm39 && OUT_DIR=. Rscript $REPO/OCM_heart/core_escape_SNPs.R
#
# These do NOT replace the SNP bed. The counter still counts every chrX SNP in
# $SNPS into x_ref/x_alt; each bed below adds ONE MORE column pair holding the
# subset of those same UMIs that fall inside (or outside) it. One pass, one set
# of UMI calls, so the partitions are exactly comparable - which is the whole
# point, since escape vs non-escape and masked vs unmasked have to be measured
# on the same molecules to mean anything.
ESCAPE_BED=$BASE/GRCm39/annotation_us_mm39_core_escape.bed
XIC_BED=$BASE/GRCm39/annotation_us_mm39_Xic500kb.bed
IMP_PAT_BED=$BASE/GRCm39/annotation_us_mm39_imprinted_pat.bed
IMP_MAT_BED=$BASE/GRCm39/annotation_us_mm39_imprinted_mat.bed
# One locus per direction on its own, as a check on the pooled sets. Both are
# concentrated - Snrpn is 82% of the paternal SNPs - so if a pooled curve falls
# short of C(d) ~ 1, only a per-locus curve separates a real per-UMI error floor
# from one contaminating locus, and getting one afterwards costs another pass.
IMP_SNRPN_BED=$BASE/GRCm39/annotation_us_mm39_imprinted_snrpn.bed
IMP_DLK1_BED=$BASE/GRCm39/annotation_us_mm39_imprinted_dlk1dio3.bed

for f in "$BAM" "$SNPS"; do
  if [ ! -f "$f" ]; then echo "Missing file: $f" >&2; exit 1; fi
done

# The subsets are optional, one at a time: a missing bed drops its columns with
# a warning rather than failing the job, so a fresh clone still produces the
# original output. Names are lower case because ase_tile_sweep.R lowercases
# column names on read, and mixed case there reads as a second set.
SUBSET_ARGS=()
add_subset() {   # name, bed, [complement]
  if [ ! -f "$2" ]; then
    echo "NOTE: no $2 - the '$1' columns will be absent." >&2
    echo "      Build the interval beds with OCM_heart/core_escape_SNPs.R." >&2
    return
  fi
  if [ "${3:-}" = "complement" ]; then
    SUBSET_ARGS+=(--subset-bed "$1=$2:complement")
  else
    SUBSET_ARGS+=(--subset-bed "$1=$2")
  fi
}
# Task 4: partition the C(d) floor into escape and non-escape chrX.
add_subset escape    "$ESCAPE_BED"
add_subset nonescape "$ESCAPE_BED" complement
# Task 2: the Xic mask, and the region it removes, side by side with the
# unmasked chrX columns rather than instead of them.
add_subset xic       "$XIC_BED"
add_subset noxic     "$XIC_BED" complement
# Task 3: the imprinted positive control. Split by imprinting direction -
# H19/Igf2, Airn/Igf2r and Meg3/Dlk1 are reciprocal pairs, and pooling them
# would put opposite alleles in one set and land C(d) near 0.5 on clean data.
# These beds name autosomes, which the counter fetches interval by interval, so
# the control costs seconds rather than a chromosome of memory.
add_subset imppat    "$IMP_PAT_BED"
add_subset impmat    "$IMP_MAT_BED"
add_subset impsnrpn  "$IMP_SNRPN_BED"
add_subset impmeg3   "$IMP_DLK1_BED"
if [ ! -d "$SPATIAL_DIR" ]; then echo "Missing dir: $SPATIAL_DIR" >&2; exit 1; fi
if [ ! -f "$BAM.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
  echo "No index next to $BAM - the counter fetches by chromosome and needs one." >&2
  exit 1
fi

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

##### ---- pass 1: per-2um-bin allelic counts (pysam) ---- #####
# -s is not enough: a failed run still leaves a header line behind, and reusing
# that silently feeds an empty table to the sweep. Require actual data rows.
#
# COMPLETENESS, not just existence. An existing table from before the subset and
# window flags were added has neither, and reusing it would run the whole
# downstream analysis on counts that cannot answer the question - silently, with
# every subset curve simply absent. So the reuse check now also demands the
# window table and every requested subset column, and re-runs the pass if any is
# missing. That does mean the first run after this change recounts both samples.
WINDOW_TSV=$OUT_DIR/chrX_windows_${WINDOW_TILE_UM}um.tsv.gz
counts_complete() {
  [ "$(head -2 "$COUNTS" 2>/dev/null | wc -l)" -gt 1 ] || return 1
  [ -s "$WINDOW_TSV" ] || { echo "  (no $WINDOW_TSV)"; return 1; }
  local hdr; hdr=$(head -1 "$COUNTS")
  local i name
  for ((i = 0; i < ${#SUBSET_ARGS[@]}; i++)); do
    [ "${SUBSET_ARGS[i]}" = "--subset-bed" ] || continue
    name=${SUBSET_ARGS[i+1]%%=*}
    case "	$hdr	" in
      *"	${name}_ref	"*) ;;
      *) echo "  (no ${name}_ref column in $COUNTS)"; return 1 ;;
    esac
  done
  return 0
}
if counts_complete && [ "$FORCE" != "1" ] && [ "$FORCE" != "force" ]; then
  echo "Reusing existing $COUNTS ($(($(wc -l < "$COUNTS") - 1)) bins; pass 'force' as arg 2 to redo)"
else
  conda activate sinto     # pysam lives here alongside sinto itself

  # Xist is excluded from the SNP set (it reports the inactive X and would
  # invert the ratio), matching the snRNA-seq runs. The autosomal control uses
  # small chromosomes on purpose: the null only needs depth comparable to chrX,
  # and peak memory scales with UMIs on the largest chromosome processed.
  # Bins come from the sb tag (sb:Z:s_002um_01679_00804-1), which spaceranger
  # writes on every read whose barcode corrected to a slide position. CB holds
  # the raw sequence and cannot be placed without the slide whitelist, and
  # outs/barcode_mappings.parquet does not carry sequences at all.
  #
  # No --require-gene: GX marks exonic assignment, while Allelome.PRO2 counts
  # over gene-body intervals from its annotation bed, introns included. On this
  # nuclear-heavy 3' data most reads are intronic (RE:A:I), so requiring GX
  # would both discard most of the data and measure a stricter thing than the
  # pipeline being sized.
  # --window-out is what makes task 1 answerable at all: the main table sums
  # over chromosome position, so nothing downstream can say WHERE on chrX the
  # excess UMIs sit. 100 kb windows on the 64um tile grid, which is the grid
  # ase_chrX_window_profile.R stratifies by X-share decile. Sparse, so the file
  # is small - only windows that carry a UMI appear.
  python3 "$REPO/spatial/ase_bin_allele_counts.py" \
    --bam "$BAM" \
    --spatial-dir "$SPATIAL_DIR" \
    --snps "$SNPS" \
    --out "$COUNTS" \
    --x-chrom chrX \
    --autosome-chroms chr19,chr18,chr17 \
    --snp-offset -1 \
    --min-mapq 255 \
    --min-baseq 20 \
    --window-out "$OUT_DIR/chrX_windows_${WINDOW_TILE_UM}um.tsv.gz" \
    --window-size 100000 \
    --window-tile-um "$WINDOW_TILE_UM" \
    --locus-out "$OUT_DIR/subset_locus_counts.tsv" \
    "${SUBSET_ARGS[@]}" \
    --threads "$SLURM_CPUS_PER_TASK"
  rc=$?
  conda deactivate
  if [ $rc -ne 0 ]; then echo "Counting pass failed" >&2; exit $rc; fi
fi

# --snp-offset -1 is not a guess: the probe scored -1 at 100.0%, +0 at 49.0%
# and +1 at 33.5% on 200k chrX reads. Column 2 of this SNP bed is a 1-based VCF
# position written as `pos, pos+1`, so it sits one to the right of the 0-based
# coordinate pysam reports. Re-probe if the bed is ever regenerated:
#   python3 $REPO/spatial/ase_bin_allele_counts.py --bam ... --out /dev/null --probe 200000

##### ---- pass 2: tile-size sweep (R) ---- #####
conda activate RNAseq
if [ -n "$TILE_UM" ]; then
  echo "Writing the sinto tile map at ${TILE_UM}um"
else
  echo "No tile size given - sweeping only, no sinto map will be written."
  echo "Re-run as: sbatch \$0 <tile_um>   (e.g. 64) once you have chosen one."
fi
SAMPLE="$SAMPLE" TILE_UM="$TILE_UM" Rscript "$REPO/spatial/ase_tile_sweep.R"
