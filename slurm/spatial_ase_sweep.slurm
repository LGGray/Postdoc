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
# exists, so tuning the sweep is just `sbatch` again. FORCE=1 redoes it.

module load slurm_setup

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
SNPS=$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.clean.bed

for f in "$BAM" "$SNPS"; do
  if [ ! -f "$f" ]; then echo "Missing file: $f" >&2; exit 1; fi
done
if [ ! -d "$SPATIAL_DIR" ]; then echo "Missing dir: $SPATIAL_DIR" >&2; exit 1; fi
if [ ! -f "$BAM.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
  echo "No index next to $BAM - the counter fetches by chromosome and needs one." >&2
  exit 1
fi

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

##### ---- pass 1: per-2um-bin allelic counts (pysam) ---- #####
if [ -s "$COUNTS" ] && [ "${FORCE:-0}" != "1" ]; then
  echo "Reusing existing $COUNTS (FORCE=1 to redo)"
else
  conda activate sinto     # pysam lives here alongside sinto itself

  # Xist is excluded from the SNP set (it reports the inactive X and would
  # invert the ratio), matching the snRNA-seq runs. The autosomal control uses
  # small chromosomes on purpose: the null only needs depth comparable to chrX,
  # and peak memory scales with UMIs on the largest chromosome processed.
  python3 "$REPO/spatial/ase_bin_allele_counts.py" \
    --bam "$BAM" \
    --spatial-dir "$SPATIAL_DIR" \
    --snps "$SNPS" \
    --out "$COUNTS" \
    --x-chrom chrX \
    --autosome-chroms chr19,chr18,chr17 \
    --min-mapq 255 \
    --min-baseq 20 \
    --threads "$SLURM_CPUS_PER_TASK"
  rc=$?
  conda deactivate
  if [ $rc -ne 0 ]; then echo "Counting pass failed" >&2; exit $rc; fi
fi

# Read the allele-match rate out of the log before trusting anything below.
# Well above 95% means the SNP bed's coordinate convention matches the BAM;
# anything near 50% is an off-by-one and --snp-offset -1 is the fix.

##### ---- pass 2: tile-size sweep (R) ---- #####
conda activate RNAseq
SAMPLE="$SAMPLE" Rscript "$REPO/spatial/ase_tile_sweep.R"
