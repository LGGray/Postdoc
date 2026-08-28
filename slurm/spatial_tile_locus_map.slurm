#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_tile_locus_map
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

# The tile ratio map without sinto and without Allelome.PRO2.
#
# WHAT THIS REPLACES
#   spatial_sinto_tiles.slurm splits the BAM into ~4000 per-tile BAMs and
#   spatial_tile_map.slurm scores each with Allelome.PRO2 - ~14 h/sample, and
#   the per-tile tables are then aggregated to whole chromosomes anyway, because
#   Allelome.PRO2 was run against chr_annotation_mm39.bed (whole-chromosome
#   intervals). ase_tile_locus_counts.py reads the BAM ONCE and accumulates per
#   (tile, chromosome) directly off the `sb` tag. Same statistic, one pass.
#
# NOTHING EXISTING IS OVERWRITTEN. Everything new lands in
#   adult_aged_spatial/ase_pysam_64um/
# except the Allelome.PRO2-format adapter tree, which has to sit where
# tile_ratio_map.R looks for it:
#   adult_aged_spatial/ase/<sample>/Allelome.PRO2_tiles_64um_pysam/
# That directory name is new, so the real Allelome.PRO2_tiles_64um tree beside
# it is untouched, as are tile_ratio_map_64um.{csv,pdf}.
#
# WHY AN ADAPTER TREE RATHER THAN A NEW PLOTTING SCRIPT
#   tile_ratio_map.R carries the lab's allelic-ratio ramp, the H&E overlay, the
#   overdispersed autosomal null and the call logic. Emitting its input format
#   means it runs UNMODIFIED, so the only difference between the two sets of
#   figures is where the counts came from - which is the comparison being made.
#   A reimplementation would have produced similar figures, not comparable ones.
#
# Array task 1 = 9w, task 2 = 78w. Whichever task finishes second draws the
# figures for both samples, so no second sbatch is needed. The lock is an
# atomic mkdir, so if the two ever finish in the same instant exactly one plots.
#
#   sbatch ~/Postdoc/slurm/spatial_tile_locus_map.slurm             # 64 um
#   sbatch ~/Postdoc/slurm/spatial_tile_locus_map.slurm 32          # 32 um
#   sbatch ~/Postdoc/slurm/spatial_tile_locus_map.slurm 64 force    # recount
#   sbatch ~/Postdoc/slurm/spatial_tile_locus_map.slurm 64 "" plot  # replot only
#
# Outputs, in adult_aged_spatial/ase_pysam_64um/:
#   <sample>/tile_chrom_counts.tsv         tile x chromosome, UMI and read level
#   <sample>/tile_gene_counts.tsv.gz       tile x gene, both levels - the thing
#                                          the Allelome.PRO2 route could not
#                                          produce per tile at all
#   <sample>/tile_locus_counts.provenance.tsv
#   tile_ratio_map_64um_pysam.pdf          the map, panels as before
#   tile_ratio_map_64um_pysam_distribution.pdf
#   tile_ratio_map_64um_pysam.csv          same schema as tile_ratio_map_64um.csv
#   tile_ratio_map_64um_pysam_provenance.tsv

module load slurm_setup

# Positional, not environment: --export=NONE means `TILE_UM=64 sbatch ...` sets
# the variable in the login shell and the job never sees it.
TILE_UM="${1:-64}"
FORCE="${2:-}"
MODE="${3:-}"           # 'plot' skips counting and only redraws

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
REPO=$BASE/Postdoc
SAMPLES="9w,78w"

# The label that names the adapter tree and the output PDFs. tile_ratio_map.R
# builds both from SNP_LABEL, so this one string keeps the new run beside the
# old one instead of on top of it.
LABEL=pysam

OUT_ROOT=$SPATIAL/ase_pysam_${TILE_UM}um
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
SNPS=$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed
OUT_DIR=$OUT_ROOT/$SAMPLE
AP2_DIR=$SPATIAL/ase/$SAMPLE/Allelome.PRO2_tiles_${TILE_UM}um_${LABEL}
COUNTS=$OUT_DIR/tile_chrom_counts.tsv

# Genome-wide gene annotation, for the optional tile x gene table. This is the
# same annotation OCM_heart/core_escape_SNPs.R builds its beds from, one row per
# transcript; Intervals.name_at() keys on the name, so transcripts of one gene
# aggregate. Optional: absent, the run still produces the tile x chromosome
# table and the map, which is what replicates Allelome.PRO2.
GENE_BED=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation_us.bed

for f in "$BAM" "$SNPS"; do
  if [ ! -f "$f" ]; then echo "Missing file: $f" >&2; exit 1; fi
done
if [ ! -d "$SPATIAL_DIR" ]; then echo "Missing dir: $SPATIAL_DIR" >&2; exit 1; fi
if [ ! -f "$BAM.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
  echo "No index next to $BAM - the counter fetches by chromosome and needs one." >&2
  exit 1
fi

GENE_ARGS=()
if [ -f "$GENE_BED" ]; then
  GENE_ARGS=(--annotation "$GENE_BED")
else
  echo "NOTE: no $GENE_BED - skipping the tile x gene table." >&2
fi

FORCE_ARGS=()
if [ "$FORCE" = "force" ] || [ "$FORCE" = "1" ]; then FORCE_ARGS=(--force); fi

mkdir -p "$OUT_DIR"

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

##### ---- pass 1: per-tile allelic counts (pysam) ---- #####
if [ "$MODE" = "plot" ]; then
  echo "MODE=plot - skipping the counting pass for $SAMPLE"
elif [ "$(head -2 "$COUNTS" 2>/dev/null | wc -l)" -gt 1 ] \
     && [ "$FORCE" != "force" ] && [ "$FORCE" != "1" ]; then
  # -s is not enough: a failed run leaves a header line behind, and reusing that
  # silently feeds an empty table to the plots.
  echo "Reusing existing $COUNTS ($(($(wc -l < "$COUNTS") - 1)) rows;"
  echo "  pass 'force' as arg 2 to recount)"
else
  conda activate sinto     # pysam lives here alongside sinto itself
  # chr1..chr19 + chrX, NOT the 3-chromosome control the bin counter uses:
  # tile_ratio_map.R pools chr1..chr19 for its autosomal null, so a narrower set
  # would change the null band rather than reproduce it.
  python3 "$REPO/spatial/ase_tile_locus_counts.py" \
    --bam "$BAM" \
    --spatial-dir "$SPATIAL_DIR" \
    --snps "$SNPS" \
    --out-dir "$OUT_DIR" \
    --ap2-dir "$AP2_DIR" \
    --tile-um "$TILE_UM" \
    --snp-offset -1 \
    --min-mapq 255 \
    --min-baseq 20 \
    --ap2-level umi \
    "${GENE_ARGS[@]}" \
    "${FORCE_ARGS[@]}" \
    --threads "$SLURM_CPUS_PER_TASK"
  rc=$?
  conda deactivate
  if [ $rc -ne 0 ]; then echo "Counting pass failed for $SAMPLE" >&2; exit $rc; fi
fi

##### ---- pass 2: the map, drawn by the existing script ---- #####
# Both samples land in one PDF, so this runs once, from whichever array task
# finishes last. mkdir is atomic on this filesystem, so exactly one task wins
# even if they finish together.
ALL_READY=1
for s in ${SAMPLES//,/ }; do
  if [ ! -s "$OUT_ROOT/$s/tile_chrom_counts.tsv" ]; then ALL_READY=0; fi
done
if [ "$ALL_READY" -ne 1 ]; then
  echo "$SAMPLE done; the other sample is still counting. It will draw the map."
  exit 0
fi
if ! mkdir "$OUT_ROOT/.plot_lock_${SLURM_ARRAY_JOB_ID:-manual}" 2>/dev/null; then
  echo "$SAMPLE done; the other task already claimed the plotting step."
  exit 0
fi

echo "Both samples counted - drawing the map from $SAMPLE's task."
conda activate seurat_env    # NOT RNAseq: tile_ratio_map.R needs arrow and png
Rscript -e 'for (p in c("arrow","png")) if (!requireNamespace(p, quietly=TRUE))
              message("note: R package \"", p, "\" not installed - a panel will be skipped")' \
  2>&1 | sed 's/^/  /'

# SCRATCH is deliberately cleared. tile_ratio_map.R supplements the collected
# tree with anything Allelome.PRO2 has left in $SCRATCH, which here would mix
# genuine read-level Allelome.PRO2 tables into the pysam figures - the one thing
# that would make the comparison meaningless. Clearing it makes that read a
# no-op; the script logs how many tables came from each source, so the zero is
# visible rather than assumed.
# COUNT_UNIT must match --ap2-level above. At 'umi' one unit in the locus tables
# is a MOLECULE, not a read, and every panel label has to say so - reads vs
# molecules is the entire difference between this map and the Allelome.PRO2 one
# (chrX pooled 0.87 per molecule against 0.68/0.60 per read on the same data).
env -u SCRATCH \
  TILE_UM="$TILE_UM" \
  SAMPLES="$SAMPLES" \
  SNP_LABEL="$LABEL" \
  COUNT_UNIT="molecule" \
  OUT_PDF="$OUT_ROOT/tile_ratio_map_${TILE_UM}um_${LABEL}.pdf" \
  Rscript "$REPO/spatial/tile_ratio_map.R"
rc=$?
if [ $rc -ne 0 ]; then echo "Plotting failed" >&2; exit $rc; fi

# The provenance sidecar tile_ratio_map.R writes looks for a SNP bed named after
# SNP_LABEL, so its snp_bed_md5 will read "bed not readable from here" - there is
# no SNPfile_..._mm39_pysam.bed. The real fingerprint is in each sample's
# tile_locus_counts.provenance.tsv, written by the counter.
echo
echo "Done. New outputs in $OUT_ROOT:"
ls -la "$OUT_ROOT"
echo
echo "Compare against the Allelome.PRO2 route:"
echo "  new: $OUT_ROOT/tile_ratio_map_${TILE_UM}um_${LABEL}.csv"
echo "  old: $SPATIAL/ase/tile_ratio_map_${TILE_UM}um.csv"
