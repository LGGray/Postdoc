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
#   sbatch ~/Postdoc/slurm/spatial_ase_sweep.slurm 64     # writes the map
#   sbatch ~/Postdoc/slurm/spatial_sinto_tiles.slurm 64   # this job
#
# Array task 1 = 9w, task 2 = 78w.

module load slurm_setup

# Positional, not environment: --export=NONE above stops `TILE_UM=64 sbatch`
# from reaching the job.
TILE_UM="${1:-${TILE_UM:-}}"
# Optional arg 2: stop after this many Allelome.PRO2 tiles. Use it to time the
# per-tile cost before committing to all ~3700 - each invocation re-parses the
# genome-wide SNP bed, and whether that is 2 minutes or 20 decides whether the
# full run fits in the wall clock.
MAX_TILES="${2:-0}"
if [ -z "$TILE_UM" ]; then
  echo "Usage: sbatch $0 <tile_um>    e.g. sbatch $0 64" >&2
  echo "Must match the size the sweep wrote a map for." >&2
  exit 1
fi

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

# NOT CB. On Visium HD the CB tag carries the raw spatial barcode sequence
# (GATGCGCTCACACTG...), which has no bin identity without the slide whitelist.
# The 2um bin is in sb (sb:Z:s_002um_01679_00804-1), and that is what the tile
# map is keyed on. Splitting on CB would match nothing and write no BAMs.
BARCODE_TAG=sb

# Genome-wide, deliberately. The autosomes are not overhead here: each tile's
# autosomal ratio is its own technical control, so a tile with a distorted
# autosomal ratio can be flagged rather than trusted. The cost is that every
# per-tile invocation parses the full SNP and annotation beds.
ANNOTATION=$BASE/GRCm39/chr_annotation_mm39.bed
# Xist reports the inactive X and would invert every ratio, so it stays out,
# matching the snRNA-seq runs in scAllelome_*.slurm.
SNPS=$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed

MIN_READ=1
TOTAL_READ=1

FINAL_DIR=$SPATIAL/ase/$SAMPLE/Allelome.PRO2_tiles_${TILE_UM}um

for f in "$BAM" "$CELLS" "$ANNOTATION" "$SNPS" "$PIPELINE"; do
  if [ ! -f "$f" ]; then echo "Missing file: $f" >&2; exit 1; fi
done
if [ -z "$SCRATCH" ]; then echo "SCRATCH is not set - aborting" >&2; exit 1; fi

mkdir -p "$FINAL_DIR"
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

# Persistent across retries, deliberately: no job id in the path. Splitting a
# ~200GB BAM into 3700 pieces is hours of work, and a job that hits the wall
# clock during Allelome.PRO2 should not have to redo it.
WORK=$SCRATCH/spatial_tiles_${SAMPLE}_${TILE_UM}um
BAM_DIR=$WORK/bam
AP2_DIR=$WORK/allelome
mkdir -p "$BAM_DIR" "$AP2_DIR"

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh

##### ---- split ---- #####
N_BAM=$(find "$BAM_DIR" -name '*.bam' | wc -l)
if [ "$N_BAM" -ge "$N_TILES" ]; then
  echo "[$(date +%T)] Reusing $N_BAM tile BAMs already in $BAM_DIR"
else
  echo "[$(date +%T)] sinto: splitting into $N_TILES tiles"
  conda activate sinto
  if ! sinto filterbarcodes --help 2>&1 | grep -q -- '--barcodetag'; then
    echo "ERROR: this sinto has no --barcodetag; it would split on CB and match" >&2
    echo "       nothing. Upgrade sinto, or pre-tag the BAM with sb copied to CB." >&2
    exit 1
  fi
  sinto filterbarcodes \
    --bam "$BAM" \
    --cells "$CELLS" \
    --outdir "$BAM_DIR" \
    --barcodetag "$BARCODE_TAG" \
    --nproc "$SLURM_CPUS_PER_TASK" || { echo "sinto failed" >&2; exit 1; }
  conda deactivate
  N_BAM=$(find "$BAM_DIR" -name '*.bam' | wc -l)
  echo "[$(date +%T)] sinto wrote $N_BAM / $N_TILES tile BAMs"
fi

if [ "$N_BAM" -eq 0 ]; then
  echo "ERROR: no BAMs written. Check that the HD BAM carries 2um bins in $BARCODE_TAG:" >&2
  echo "       samtools view $BAM | head -1 | tr '\\t' '\\n' | grep '^$BARCODE_TAG:'" >&2
  exit 1
fi

##### ---- Allelome.PRO2, one invocation per tile ---- #####
conda activate RNAseq
module load parallel

# Only tiles without a locus table. A timed-out run therefore picks up where it
# stopped instead of starting over.
: > "$WORK/task-list.txt"
n_todo=0; n_done=0
while IFS= read -r bam; do
  tile=$(basename "$bam" .bam)
  if [ -s "$AP2_DIR/$tile/locus_table.txt" ] || [ -s "$FINAL_DIR/$tile/locus_table.txt" ]; then
    n_done=$((n_done + 1)); continue
  fi
  echo "$PIPELINE -i $bam -a $ANNOTATION -s $SNPS -r $MIN_READ -t $TOTAL_READ -o $AP2_DIR"
  n_todo=$((n_todo + 1))
done < <(find "$BAM_DIR" -name '*.bam' | sort) >> "$WORK/task-list.txt"

if [ "$MAX_TILES" -gt 0 ]; then
  head -n "$MAX_TILES" "$WORK/task-list.txt" > "$WORK/task-list.head.txt"
  mv "$WORK/task-list.head.txt" "$WORK/task-list.txt"
  n_todo=$(wc -l < "$WORK/task-list.txt")
  echo "TIMING RUN: limited to $n_todo tiles"
fi
echo "[$(date +%T)] Allelome.PRO2: $n_todo to run, $n_done already done"

nproc=$((SLURM_CPUS_PER_TASK - 2))
if [ "$nproc" -lt 1 ]; then nproc=1; fi
start=$SECONDS
parallel -j "$nproc" --joblog "$WORK/parallel.log" :::: "$WORK/task-list.txt"
elapsed=$((SECONDS - start))

if [ "$n_todo" -gt 0 ]; then
  per_tile=$(awk -v e="$elapsed" -v n="$n_todo" -v j="$nproc" \
                 'BEGIN{printf "%.1f", e*j/n}')
  total_h=$(awk -v p="$per_tile" -v n="$N_TILES" -v j="$nproc" \
                'BEGIN{printf "%.1f", p*n/j/3600}')
  echo "[$(date +%T)] $n_todo tiles in ${elapsed}s at -j $nproc"
  echo "  ~${per_tile}s of CPU per tile -> ~${total_h}h for all $N_TILES tiles"
fi

FAILED=$(awk 'NR>1 && $7 != 0' "$WORK/parallel.log" 2>/dev/null | wc -l)
echo "Allelome.PRO2 failures this pass: $FAILED / $n_todo"

##### ---- collect ---- #####
# Only the locus tables and BED files are worth keeping; the tile BAMs are
# reproducible from the map and would be tens of thousands of small files on
# permanent storage.
mkdir -p "$FINAL_DIR"
shopt -s nullglob
for d in "$AP2_DIR"/*/; do
  name=$(basename "$d")
  mkdir -p "$FINAL_DIR/$name"
  [ -f "$d/locus_table.txt" ] && mv "$d/locus_table.txt" "$FINAL_DIR/$name/"
  [ -d "$d/BED_files" ] && mv "$d"/BED_files/*.bed "$FINAL_DIR/$name/" 2>/dev/null
done
cp "$WORK/parallel.log" "$FINAL_DIR/" 2>/dev/null

n_final=$(find "$FINAL_DIR" -name locus_table.txt | wc -l)
echo "Results in $FINAL_DIR: $n_final / $N_TILES tiles have a locus table"
if [ "$n_final" -lt "$N_TILES" ] && [ "$MAX_TILES" -eq 0 ]; then
  echo "Incomplete - resubmit the same command and it will run only what is left."
fi
echo "Tile BAMs kept in $BAM_DIR - delete once you are happy with the output."
