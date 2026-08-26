#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_sinto_tiles
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --array=1-2

# Split the Visium HD BAM into square tiles and run Allelome.PRO2 on each.
#
# Run this only after spatial_ase_sweep.slurm and only at the tile size that
# sweep chose - the barcode -> tile map it consumes is written by the sweep's R
# step when TILE_UM is set:
#
#   TILE_UM=64 sbatch ~/Postdoc/slurm/spatial_ase_sweep.slurm     # writes the map
#   TILE_UM=64 sbatch ~/Postdoc/slurm/spatial_sinto_tiles.slurm   # this job
#
# Array task 1 = 9w, task 2 = 78w.

module load slurm_setup

TILE_UM=${TILE_UM:?Set TILE_UM to the size chosen by the sweep, e.g. TILE_UM=64 sbatch ...}

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
PIPELINE=$BASE/Allelome.LINK/00_src/Allelome.PRO2.sh

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
  SAMPLE=9w
elif [ "$SLURM_ARRAY_TASK_ID" -eq 2 ]; then
  SAMPLE=78w
else
  echo "Unexpected SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2; exit 1
fi

BAM=$SPATIAL/$SAMPLE/outs/possorted_genome_bam.bam
CELLS=$SPATIAL/ase/$SAMPLE/sinto_tiles_${TILE_UM}um.tsv

# chrX-only annotation keeps Allelome.PRO2 fast when it is being invoked once
# per tile. Swap to chr_annotation_mm39.bed for a genome-wide run.
ANNOTATION=$BASE/GRCm39/chr_annotation_mm39.bed
# Xist reports the inactive X and would invert every ratio, so it stays out,
# matching the snRNA-seq runs in scAllelome_*.slurm.
SNPS=$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_withoutXist.clean.bed

MIN_READ=1
TOTAL_READ=1

FINAL_DIR=$SPATIAL/ase/$SAMPLE/Allelome.PRO2_tiles_${TILE_UM}um

for f in "$BAM" "$CELLS" "$ANNOTATION" "$SNPS" "$PIPELINE"; do
  if [ ! -f "$f" ]; then echo "Missing file: $f" >&2; exit 1; fi
done
if [ -z "$SCRATCH" ]; then echo "SCRATCH is not set - aborting" >&2; exit 1; fi

N_TILES=$(cut -f2 "$CELLS" | sort -u | wc -l)
echo "Tile size ${TILE_UM}um, $N_TILES tiles, $(wc -l < "$CELLS") bins"

# sinto filterbarcodes keeps one output BAM open per group, so the tile count
# runs straight into the open-file limit. Raise it as far as the hard limit
# allows and refuse to start if that is still not enough - failing here is much
# cheaper than discovering it eight hours in with half the tiles truncated.
HARD=$(ulimit -Hn)
ulimit -n "$HARD" 2>/dev/null || true
echo "open file limit: $(ulimit -n) (hard $HARD)"
if [ "$N_TILES" -gt $(( $(ulimit -n) - 256 )) ]; then
  echo "ERROR: $N_TILES tiles exceeds the open-file budget. Raise TILE_UM or" >&2
  echo "       SINTO_MIN_UMI in ase_tile_sweep.R and regenerate the map." >&2
  exit 1
fi

SCRATCH_DIR=$SCRATCH/spatial_tiles_${SAMPLE}_${TILE_UM}um_${SLURM_JOB_ID}
BAM_DIR=$SCRATCH_DIR/bam
AP2_DIR=$SCRATCH_DIR/allelome
rm -rf "$SCRATCH_DIR"
mkdir -p "$BAM_DIR" "$AP2_DIR"

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

##### ---- split ---- #####
conda activate sinto
sinto filterbarcodes \
  --bam "$BAM" \
  --cells "$CELLS" \
  --outdir "$BAM_DIR" \
  --barcodetag CB \
  --nproc "$SLURM_CPUS_PER_TASK" || { echo "sinto failed" >&2; exit 1; }
conda deactivate

N_BAM=$(find "$BAM_DIR" -name '*.bam' | wc -l)
echo "sinto wrote $N_BAM / $N_TILES tile BAMs"
if [ "$N_BAM" -eq 0 ]; then
  echo "ERROR: no BAMs written. Check that the HD BAM carries 2um barcodes in CB:" >&2
  echo "       samtools view $BAM | head -1 | tr '\\t' '\\n' | grep '^CB:'" >&2
  exit 1
fi

##### ---- Allelome.PRO2, one invocation per tile ---- #####
conda activate RNAseq
module load parallel

find "$BAM_DIR" -name '*.bam' | sort > "$SCRATCH_DIR/bamfiles.txt"
while IFS= read -r bam; do
  echo "$PIPELINE -i $bam -a $ANNOTATION -s $SNPS -r $MIN_READ -t $TOTAL_READ -o $AP2_DIR"
done < "$SCRATCH_DIR/bamfiles.txt" > "$SCRATCH_DIR/task-list.txt"

nproc=$((SLURM_CPUS_PER_TASK - 2))
if [ "$nproc" -lt 1 ]; then nproc=1; fi
parallel -j "$nproc" --joblog "$SCRATCH_DIR/parallel.log" :::: "$SCRATCH_DIR/task-list.txt"

FAILED=$(awk 'NR>1 && $7 != 0' "$SCRATCH_DIR/parallel.log" | wc -l)
echo "Allelome.PRO2: $FAILED / $(wc -l < "$SCRATCH_DIR/task-list.txt") tiles failed"

##### ---- collect ---- #####
# Only the locus tables and BED files are worth keeping; the tile BAMs are
# reproducible from the map and would be tens of thousands of small files on
# permanent storage.
mkdir -p "$FINAL_DIR"
for d in "$AP2_DIR"/*/; do
  name=$(basename "$d")
  mkdir -p "$FINAL_DIR/$name"
  [ -f "$d/locus_table.txt" ] && mv "$d/locus_table.txt" "$FINAL_DIR/$name/"
  [ -d "$d/BED_files" ] && mv "$d"/BED_files/*.bed "$FINAL_DIR/$name/" 2>/dev/null
done
cp "$SCRATCH_DIR/parallel.log" "$FINAL_DIR/" 2>/dev/null

echo "Results in $FINAL_DIR"
echo "Tile BAMs left in $BAM_DIR - delete once you are happy with the output."
