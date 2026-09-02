#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_tile_nested_loci
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
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --array=1-2

# Recount ONE chromosome so a nested gene can be seen at all.
#
# THE PROBLEM THIS SOLVES. Intervals.name_at() in ase_bin_allele_counts.py gives
# a SNP to the containing interval with the SMALLEST START - one gene per
# position, deterministically. Where one gene sits inside another, the outer one
# takes every SNP and the inner one reads zero molecules everywhere. Kcnq1ot1
# (chr7:142,766,847-142,850,284) is entirely inside Kcnq1
# (chr7:142,660,613-142,980,787), which starts 106 kb earlier, so the main tree
# reports 905 molecules for Kcnq1 and none for Kcnq1ot1.
#
# That is not only a missing gene. Kcnq1 is maternally expressed and Kcnq1ot1 is
# its PATERNALLY expressed antisense over the same coordinates, so what the main
# tree calls "Kcnq1" is the reciprocal pair summed - AR 0.602 at 9w and 0.494 at
# 78w, where a clean maternal control reads ~1. The imprinted control panel was
# showing an average of two opposite directions.
#
# WHAT THIS RUNS. The same counter, same conventions, same SNP mask, restricted
# to chr7 and pointed at a bed holding ONLY the nested genes' own intervals. On
# that bed no earlier interval exists to take the SNPs, so they are attributed
# where they belong. Everything lands in a SIDE TREE:
#     adult_aged_spatial/ase_pysam_dup_64um_nested/<sample>/
# and the figures merge it with EXTRA_ROOT, which REPLACES the same-named rows
# from the main tree rather than adding to them. The main tree is not touched.
#
# COST. chr7 alone, not the genome: the full 20-chromosome pass took ~4 h per
# sample, so this is well inside the 6 h wall time with room for a slow node.
# Both samples run as one array.
#
#   sbatch ~/Postdoc/slurm/spatial_tile_nested_loci.slurm
#   sbatch ~/Postdoc/slurm/spatial_tile_nested_loci.slurm 64 force   # recount
#
# READ THE COUNTS BEFORE THE FIGURES. Kcnq1ot1 is paternally expressed, so with
# B6 mother x CAST father it should read AR ~0 (CAST). If it comes back near
# 0.5, the SNPs are not being re-attributed and the bed is wrong - check the
# gene_rows line in the provenance file before drawing anything.

module load slurm_setup

TILE_UM="${1:-64}"
FORCE="${2:-}"

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
REPO=$BASE/Postdoc
REF=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39
ANNOT=$REF/20250512_RefSeq/annotation_us.bed
SNPS=$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed

# The nested genes to rescue, and the chromosome they live on. One chromosome
# per run keeps this cheap; add a gene here only if it is on CHROM.
GENES="Kcnq1ot1"
CHROM=chr7

OUT_ROOT=$SPATIAL/ase_pysam_dup_${TILE_UM}um_nested
BED=$OUT_ROOT/annotation_us_mm39_nested_${CHROM}.bed
mkdir -p "$OUT_ROOT"

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
  SAMPLE=9w
elif [ "$SLURM_ARRAY_TASK_ID" -eq 2 ]; then
  SAMPLE=78w
else
  echo "Unexpected SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2; exit 1
fi

BAM=$SPATIAL/$SAMPLE/outs/possorted_genome_bam.bam
SPATIAL_DIR=$SPATIAL/$SAMPLE/outs/binned_outputs/square_002um/spatial
OUT_DIR=$OUT_ROOT/$SAMPLE

for f in "$BAM" "$SNPS" "$ANNOT"; do
  if [ ! -f "$f" ]; then echo "Missing file: $f" >&2; exit 1; fi
done
if [ ! -d "$SPATIAL_DIR" ]; then echo "Missing dir: $SPATIAL_DIR" >&2; exit 1; fi

# The bed: one row per transcript of the named genes, nothing else. Built here
# rather than checked in, because it has to match whichever annotation_us.bed
# the main tree was counted against - a stale copy would re-attribute the SNPs
# to slightly different spans and the merge would be comparing two annotations.
# Task 1 writes it; task 2 waits, so the two do not race on the same path.
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ] || [ ! -s "$BED" ]; then
  awk -v genes="$GENES" -v chrom="$CHROM" '
    BEGIN { n = split(genes, g, ","); for (i = 1; i <= n; i++) want[g[i]] = 1 }
    $1 == chrom && ($4 in want) { print; seen[$4] = 1 }
    END { for (k in want) if (!(k in seen))
            printf "NOT FOUND on %s: %s\n", chrom, k > "/dev/stderr" }
  ' "$ANNOT" > "$BED.tmp" && mv "$BED.tmp" "$BED"
fi
for _ in $(seq 1 60); do [ -s "$BED" ] && break; sleep 5; done
if [ ! -s "$BED" ]; then
  echo "Empty bed $BED - none of '$GENES' is on $CHROM in $ANNOT" >&2; exit 1
fi
echo "bed: $BED"
awk '{printf "  %s:%d-%d  %s\n", $1, $2, $3, $4}' "$BED"
echo "  md5 $(md5sum "$BED" | cut -d" " -f1)"

FORCE_ARGS=()
if [ "$FORCE" = "force" ] || [ "$FORCE" = "1" ]; then FORCE_ARGS=(--force); fi

mkdir -p "$OUT_DIR"
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate sinto      # pysam lives here alongside sinto itself

# Every flag matches spatial_tile_locus_map.slurm's dup run except --chroms and
# --annotation. They have to: the whole point is that these counts drop into the
# same figures as rows of the same kind, and a different mapq or offset would
# make them a different measurement wearing the same column names.
python3 "$REPO/spatial/ase_tile_locus_counts.py" \
  --bam "$BAM" \
  --spatial-dir "$SPATIAL_DIR" \
  --snps "$SNPS" \
  --out-dir "$OUT_DIR" \
  --tile-um "$TILE_UM" \
  --chroms "$CHROM" \
  --annotation "$BED" \
  --snp-offset -1 \
  --min-mapq 255 \
  --min-baseq 20 \
  --ap2-level reads \
  --keep-duplicates \
  "${FORCE_ARGS[@]}" \
  --threads "$SLURM_CPUS_PER_TASK"
rc=$?
conda deactivate
if [ $rc -ne 0 ]; then echo "ase_tile_locus_counts.py failed (rc=$rc)" >&2; exit $rc; fi

echo
echo "=== what came out, per gene ==="
# Printed here rather than left for the figures: if the re-attribution did not
# work, this is where it is cheap to see. Kcnq1ot1 is paternal, so a1 (B6)
# should be the SMALL side.
gzip -cd "$OUT_DIR/tile_gene_counts.tsv.gz" 2>/dev/null | awk -F'\t' '
  NR == 1 { next }
  { n[$5] += $8; a1[$5] += $6; t[$5]++ }
  END { printf "  %-12s %8s %8s %7s\n", "gene", "umi", "tiles", "AR"
        for (g in n) printf "  %-12s %8d %8d %7.3f\n", g, n[g], t[g], a1[g]/n[g] }'

echo
echo "Done. Merge into the figures with:"
echo "  EXTRA_ROOT=$OUT_ROOT"
