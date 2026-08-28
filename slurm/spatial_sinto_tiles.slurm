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
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --array=1-2

# 12h, sized from measurement rather than guessed. Job 207291 did 701/3738
# tiles (9w) and 711/3215 (78w) in 2h35 at -j 30, i.e. ~270 tiles/h and 6.6 min
# of CPU per tile - that cost is the per-invocation genome-wide SNP and
# annotation parse, so it is fixed per tile and does not fall as the run
# proceeds. Whole-sample totals are ~414 CPU-h (9w) and ~351 (78w), which is
# 13.8h and 11.7h of wall clock on one 32-core node. 12h therefore clears what
# is LEFT after 207291 with margin, but a cold start on 9w needs ~14h - raise
# this to 16h if the scratch tile BAMs are ever rebuilt from nothing.
#
# An earlier note here claimed the wall clock was the queue lever and cut this
# to 4h on that basis. Lachlan's read is that the memory request was the actual
# blocker, and 207291 does support that: it started promptly at 120G. The 24h
# cm4_tiny ceiling remains the hard cap if 12h ever proves too short.
#
# Nothing is lost by an under-length job in any case - every stage resumes, so a
# kill costs one resubmission and no recomputation. That guarantee is only as
# good as the resume check, which is why the locus tables are now validated for
# completeness rather than mere existence; see MIN_LOCI below.
#
# 120G rather than 150G because 32 cores on a cm4 node carry a proportional
# share of ~145G, and a request above that needs a node with more free memory
# than free cores - a much smaller target than any node with 32 idle cores.
# 48G was tried and job 207231 died on SIGKILL: 30 concurrent Allelome.PRO2
# each parse the genome-wide SNP bed, so the footprint scales with -j, not with
# the tile size, and the single-invocation 8G in Allelome.PRO.slurm does not
# bound it. Do not trim this again without a MaxRSS to point at:
#   sacct -j <jobid> -o JobID,ReqMem,MaxRSS,Elapsed,State,ExitCode
#
# The tile BAMs are already split and kept in $SCRATCH, so the defaults above
# are all that is needed from here. If the split ever has to be redone, it
# wants the opposite shape - it is I/O bound and sinto holds one output BAM
# handle per tile in every worker, so cores buy little and cost memory:
#
#   sbatch --cpus-per-task=8 --mem=100G --time=08:00:00 \
#          ~/Postdoc/slurm/spatial_sinto_tiles.slurm 64

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
# Optional arg 3: which SNP mask to score against. The bed is
# SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_<label>.bed, and every output path
# carries the label so two masks coexist instead of overwriting each other.
#
# The expensive stages - the sinto split and the duplicate/MAPQ filter - are
# MASK-INDEPENDENT and stay shared, so a second mask costs only the
# Allelome.PRO2 pass, not another split of a 200 GB BAM.
#
#   sbatch $0 64                    # no_Xist, the original paths, unchanged
#   sbatch $0 64 0 no_Xic500kb      # Xist +- 500 kb masked, alongside it
#
# Do this only if spatial_mask_comparison.slurm says the Xic carries the effect.
# It is ~14h per sample; the per-allele rates answer the same question for free.
SNP_LABEL="${3:-no_Xist}"
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
# matching the snRNA-seq runs in scAllelome_*.slurm. no_Xic500kb additionally
# drops Ftx, Jpx, Tsix and the intergenic Xic, which the Xist gene span leaves
# behind and which are effectively CAST-only in this design.
SNPS=$BASE/GRCm39/SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_${SNP_LABEL}.bed
# Empty for the default mask, so every existing path is byte-identical.
if [ "$SNP_LABEL" = "no_Xist" ]; then SUF=""; else SUF="_$SNP_LABEL"; fi

MIN_READ=1
TOTAL_READ=1

# Drop PCR duplicates and multimappers from each tile BAM before scoring.
#
# Allelome.PRO2 counts reads, and ~54% of this BAM is flagged duplicate. Every
# read of one UMI carries the same allele, so duplicates add no evidence about
# which allele a molecule came from - they only add weight. A 64um tile holds
# ~14 independent chrX molecules, so one heavily amplified molecule can move the
# tile's ratio by 0.7. Measured on 78w tiles: read-level ratios of 0.09-0.20
# where the same tiles were 0.62-0.82 at molecule level.
#
# -F 0x400 = not a duplicate; -q 255 = STAR-unique. Both match the filters the
# independent counter used, so the two methods become comparable.
DEDUP=1
FILTER_MAPQ=255

FINAL_DIR=$SPATIAL/ase/$SAMPLE/Allelome.PRO2_tiles_${TILE_UM}um${SUF}

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
# WORK and BAM_DIR carry no mask label: the tile BAMs are the same whichever
# SNP bed scores them, and re-splitting for a second mask would be hours wasted.
WORK=$SCRATCH/spatial_tiles_${SAMPLE}_${TILE_UM}um
BAM_DIR=$WORK/bam
AP2_DIR=$WORK/allelome${SUF}
mkdir -p "$BAM_DIR" "$AP2_DIR"
echo "SNP mask: $SNP_LABEL"
echo "  bed        $SNPS"
echo "  scored to  $AP2_DIR"
echo "  collected  $FINAL_DIR"

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

conda activate RNAseq
module load parallel
nproc=$((SLURM_CPUS_PER_TASK - 2))
if [ "$nproc" -lt 1 ]; then nproc=1; fi

##### ---- filter duplicates and multimappers ---- #####
if [ "$DEDUP" = "1" ]; then
  FLT_DIR=$WORK/bam_flt
  mkdir -p "$FLT_DIR"

  # Discard anything the previous pass left half-written before deciding what
  # still needs filtering. Job 207231 was killed with SIGKILL mid-stage, and a
  # BAM that was being written at that moment is non-empty but truncated - the
  # old `[ -s ]` test called it done, so Allelome.PRO2 would have scored a
  # partial tile and reported a plausible but wrong ratio, silently. quickcheck
  # reads the header and the BGZF EOF block, which is exactly what a truncated
  # file is missing. One invocation over all tiles: it is header-only, so this
  # costs seconds where 3700 separate calls would not.
  BAD=$WORK/flt-truncated.txt
  : > "$BAD"
  if compgen -G "$FLT_DIR/*.bam" > /dev/null; then
    samtools quickcheck -v "$FLT_DIR"/*.bam > "$BAD" 2>/dev/null || true
    n_bad=$(wc -l < "$BAD")
    if [ "$n_bad" -gt 0 ]; then
      echo "[$(date +%T)] discarding $n_bad truncated tile BAMs from the previous pass"
      while IFS= read -r b; do rm -f "$b" "$b.bai"; done < "$BAD"
    fi
  fi

  : > "$WORK/filter-list.txt"
  for bam in "$BAM_DIR"/*.bam; do
    t=$(basename "$bam")
    # Both, not just the BAM: the index is only written once samtools view has
    # closed the output, so its presence is a second, independent witness that
    # the pair is complete.
    [ -s "$FLT_DIR/$t" ] && [ -s "$FLT_DIR/$t.bai" ] && continue
    echo "samtools view -b -F 0x400 -q $FILTER_MAPQ -o '$FLT_DIR/$t' '$bam' && samtools index '$FLT_DIR/$t'"
  done >> "$WORK/filter-list.txt"
  n_flt=$(wc -l < "$WORK/filter-list.txt")
  if [ "$n_flt" -gt 0 ]; then
    echo "[$(date +%T)] filtering $n_flt tile BAMs (duplicates, MAPQ < $FILTER_MAPQ)"
    parallel -j "$nproc" :::: "$WORK/filter-list.txt"
  else
    echo "[$(date +%T)] filtered tile BAMs already present"
  fi
  SCORE_DIR=$FLT_DIR
  # Worth eyeballing once: this should be roughly half, and the tile ratios
  # should stop swinging between 0.1 and 1.0 once it is applied.
  ex=$(ls "$BAM_DIR"/*.bam | head -1); exb=$(basename "$ex")
  echo "  example $exb: $(samtools view -c "$ex") reads -> $(samtools view -c "$FLT_DIR/$exb") after filtering"
else
  SCORE_DIR=$BAM_DIR
  echo "WARNING: DEDUP=0, scoring raw tile BAMs. Per-tile ratios will be driven"
  echo "         by PCR duplication, not by allele composition."
fi

##### ---- Allelome.PRO2, one invocation per tile ---- #####

# Allelome.PRO2 does not write to <outdir>/<sample>/. It creates a dated job
# directory per run:
#
#   <outdir>/2026_08_26_tile_64um_r0015_c0040_chr_annotation_mm39.bed_1/
#            ^date       ^bam basename        ^annotation basename    ^run no.
#
# The trailing counter increments on every re-run, so without a resume check
# a resubmitted job silently accumulates _2, _3, ... duplicates of work that
# is already done.
ANNOT_BASE=$(basename "$ANNOTATION")

# A locus table counts as done only if it is STRUCTURALLY COMPLETE.
#
# The old test was `[ -s ]` on the collected copy, plus a bare `compgen -G` on
# the scratch copy that did not even check for non-empty. Both pass for a table
# half-written by a worker killed at the wall clock: the next run skips the tile
# and reports the truncated ratio as real, silently. Up to $nproc tiles are in
# flight at every kill, so across several resubmissions that is tens of wrong
# tiles, invisible in the maps. This is the same failure already guarded for the
# tile BAMs with `samtools quickcheck` - the locus tables never got an
# equivalent.
#
# Three checks: a header naming the columns the R side reads; every row carrying
# the header's field count, which is what a mid-line truncation breaks; and at
# least MIN_LOCI rows, since a genome-wide annotation over a tile with thousands
# of autosomal reads yields about one row per chromosome, so a table with a
# handful was cut short at a line boundary where the field count still matches.
MIN_LOCI=10
# How many times a tile carrying no chrX row may be re-run before its table is
# accepted as-is. A missing chrX row is almost always a tail truncation, but not
# necessarily: a tile can legitimately end up with no chrX read after
# -F 0x400 -q 255 even though the sweep counted >= SINTO_MIN_UMI chrX molecules,
# because that counter works on molecules and does not apply the MAPQ filter.
# Without a budget such a tile would be purged and re-run on every resubmission
# for the life of the project, at 6.6 min of CPU each time, and the tile count
# would never reach N_TILES.
MAX_REDO=2
REDO_DIR=$WORK/redo-count${SUF}
mkdir -p "$REDO_DIR"

# Exit 0 = complete. 2 = structurally sound but carrying no chrX row.
# 1 = structurally broken.
#
# chrX is the sharpest of these checks, not an extra one. chr_annotation_mm39.bed
# is sorted, so chrX follows every autosome and is the FIRST row a cut tail
# loses: a table truncated from 21 rows to 15 still has a clean field count and
# still clears MIN_LOCI, and chrX is precisely what is gone. MIN_LOCI stays as
# the floor for the opposite case, a table cut down to almost nothing.
locus_table_state() {
  [ -s "$1" ] || return 1
  awk -v min="$MIN_LOCI" '
    NR == 1 { nf = NF; if (tolower($0) !~ /total_reads/) { bad = 1; exit } next }
    NF != nf { bad = 1; exit }
    $1 == "chrX" { seen_x = 1 }
    END {
      if (bad || NR - 1 < min) exit 1
      if (!seen_x) exit 2
      exit 0
    }
  ' "$1" 2>/dev/null
}

# The tile name, from either layout:
#   $FINAL_DIR/<tile>/locus_table.txt
#   $AP2_DIR/<YYYY_MM_DD>_<tile>_<annotation>_<n>/locus_table.txt
tile_of_lt() {
  local d
  d=$(basename "$(dirname "$1")")
  d=${d#????_??_??_}
  d=${d%_${ANNOT_BASE}_*}
  printf '%s' "$d"
}

# Side-effect free, so it is safe in the resume check, the collect stage and the
# final tally. A chrX-free table is acceptable only once its redo budget is
# spent - the purge pass below is what spends it.
locus_table_ok() {
  local st tile n=0
  locus_table_state "$1"; st=$?
  [ "$st" -eq 0 ] && return 0
  [ "$st" -eq 1 ] && return 1
  tile=$(tile_of_lt "$1")
  [ -f "$REDO_DIR/$tile" ] && n=$(cat "$REDO_DIR/$tile")
  [ "$n" -ge "$MAX_REDO" ] && return 0
  return 1
}

# Copy validated locus tables and BED files from $AP2_DIR into $FINAL_DIR.
#
# Called TWICE: once here at startup and once after the parallel run. The
# startup call is the important one. This job is killed at the wall clock -
# 207291 was - and a kill takes `parallel` with it, so control never reaches the
# end-of-job call and NOTHING is copied. That is why $FINAL_DIR was empty while
# ~1400 locus tables sat in $SCRATCH: not one job had ever finished normally.
# $SCRATCH is purged on an access-time policy, so those results were only ever
# one purge away from being lost. Collecting at startup means each resubmission
# rescues the previous kill's output, and the only work ever exposed is the
# final job's - which, if it exits normally, collects itself.
#
# Only the locus tables and BED files are worth keeping; the tile BAMs are
# reproducible from the map and would be thousands of small files on permanent
# storage. The tile name has to be recovered from the dated job directory:
#   2026_08_26_tile_64um_r0015_c0040_chr_annotation_mm39.bed_1
#   -> strip the leading YYYY_MM_DD_, then the trailing _<annotation>_<run no.>
collect_tables() {
  local label="${1:-}"
  # nullglob locally: the loops below want an empty iteration rather than an
  # unmatched literal pattern, but the rest of the script has not been audited
  # under it, so restore whatever was set.
  local had_nullglob=0
  shopt -q nullglob && had_nullglob=1
  shopt -s nullglob

  local n_moved=0 n_bad_lt=0 n_already=0 d name tile lt b
  for d in "$AP2_DIR"/*_"$ANNOT_BASE"_*/; do
    name=$(basename "$d")
    tile=${name#????_??_??_}
    tile=${tile%_${ANNOT_BASE}_*}
    # Already collected and complete: skip. Without this the startup call would
    # re-copy every tile collected by every previous job, every time.
    if locus_table_ok "$FINAL_DIR/$tile/locus_table.txt"; then
      n_already=$((n_already + 1)); continue
    fi
    lt=$(find "$d" -maxdepth 2 -name locus_table.txt | head -1)
    [ -z "$lt" ] && continue
    # Never promote an incomplete table to permanent storage. Without this, a
    # table truncated after the purge pass - i.e. written by a worker killed
    # during THIS job - would be copied to $FINAL_DIR and pass every later check.
    if ! locus_table_ok "$lt"; then
      n_bad_lt=$((n_bad_lt + 1)); continue
    fi
    mkdir -p "$FINAL_DIR/$tile"
    # Write-then-rename, because cp is not atomic: a kill during the copy leaves
    # a partial file in $FINAL_DIR under the name the resume check trusts.
    # rename within one filesystem is atomic, so the destination is either the
    # old contents or the complete new ones, never half of either.
    cp "$lt" "$FINAL_DIR/$tile/.locus_table.txt.tmp" &&
      mv -f "$FINAL_DIR/$tile/.locus_table.txt.tmp" "$FINAL_DIR/$tile/locus_table.txt"
    for b in "$d"BED_files/*.bed; do
      cp "$b" "$FINAL_DIR/$tile/.$(basename "$b").tmp" &&
        mv -f "$FINAL_DIR/$tile/.$(basename "$b").tmp" "$FINAL_DIR/$tile/$(basename "$b")"
    done
    n_moved=$((n_moved + 1))
  done

  echo "[$(date +%T)] collect${label:+ ($label)}: $n_moved newly copied to $FINAL_DIR, $n_already already there"
  if [ "$n_bad_lt" -gt 0 ]; then
    echo "  skipped $n_bad_lt incomplete locus tables - those tiles stay on the to-do list"
  fi
  [ "$had_nullglob" -eq 1 ] || shopt -u nullglob
}

# Rescue whatever the previous job produced before it was killed, BEFORE
# building the to-do list - so the resume check below finds those tiles in
# $FINAL_DIR and does not re-run them.
collect_tables startup

# Discard anything a previous kill left half-written, in both locations, before
# deciding what still needs running. Removing rather than ignoring matters for
# the scratch copies: the collect stage globs that directory, so a truncated
# table left in place could still be copied over a good one.
DISCARDED=$WORK/locus-tables-discarded${SUF}.txt
n_trunc=0; n_nox=0; n_nox_kept=0
for lt in "$FINAL_DIR"/*/locus_table.txt \
          "$AP2_DIR"/*_"$ANNOT_BASE"_*/locus_table.txt; do
  [ -e "$lt" ] || continue
  locus_table_state "$lt"; st=$?
  [ "$st" -eq 0 ] && continue
  tile=$(tile_of_lt "$lt")
  if [ "$st" -eq 2 ]; then
    n=0; [ -f "$REDO_DIR/$tile" ] && n=$(cat "$REDO_DIR/$tile")
    if [ "$n" -ge "$MAX_REDO" ]; then
      # Re-run MAX_REDO times and still no chrX. Take it at face value: this
      # tile has no informative chrX read that survives the filters.
      echo "$(date +%F_%T) no-chrX-accepted-after-$n-redos  $lt" >> "$DISCARDED"
      n_nox_kept=$((n_nox_kept + 1)); continue
    fi
    echo $((n + 1)) > "$REDO_DIR/$tile"
    echo "$(date +%F_%T) no-chrX redo-$((n + 1)) $(wc -l < "$lt") lines  $lt" >> "$DISCARDED"
    n_nox=$((n_nox + 1))
  else
    echo "$(date +%F_%T) truncated $(wc -l < "$lt") lines  $lt" >> "$DISCARDED"
    n_trunc=$((n_trunc + 1))
  fi
  rm -f "$lt"
done
if [ $((n_trunc + n_nox + n_nox_kept)) -gt 0 ]; then
  echo "[$(date +%T)] locus tables from a previous pass:"
  echo "  $n_trunc structurally truncated - discarded, will be re-run"
  echo "  $n_nox complete but with no chrX row - discarded, will be re-run"
  echo "  $n_nox_kept with no chrX row after $MAX_REDO redos - accepted as genuinely chrX-free"
  echo "  details in $DISCARDED"
fi

# Only tiles without a locus table. A timed-out run therefore picks up where it
# stopped instead of starting over.
TASK_LIST=$WORK/task-list${SUF}.txt
: > "$TASK_LIST"
n_todo=0; n_done=0
while IFS= read -r bam; do
  tile=$(basename "$bam" .bam)
  if locus_table_ok "$FINAL_DIR/$tile/locus_table.txt"; then
    n_done=$((n_done + 1)); continue
  fi
  # A resubmitted tile leaves _1.._n dated directories; any one complete table
  # is enough to call the tile done.
  tile_done=0
  for lt in "$AP2_DIR"/*_"${tile}"_"${ANNOT_BASE}"_*/locus_table.txt; do
    [ -e "$lt" ] || continue
    if locus_table_ok "$lt"; then tile_done=1; break; fi
  done
  if [ "$tile_done" -eq 1 ]; then
    n_done=$((n_done + 1)); continue
  fi
  echo "$PIPELINE -i $bam -a $ANNOTATION -s $SNPS -r $MIN_READ -t $TOTAL_READ -o $AP2_DIR"
  n_todo=$((n_todo + 1))
done < <(find "$SCORE_DIR" -name '*.bam' | sort) >> "$TASK_LIST"

if [ "$MAX_TILES" -gt 0 ]; then
  head -n "$MAX_TILES" "$TASK_LIST" > "$TASK_LIST.head"
  mv "$TASK_LIST.head" "$TASK_LIST"
  n_todo=$(wc -l < "$TASK_LIST")
  echo "TIMING RUN: limited to $n_todo tiles"
fi
echo "[$(date +%T)] Allelome.PRO2: $n_todo to run, $n_done already done"

start=$SECONDS
parallel -j "$nproc" --joblog "$WORK/parallel${SUF}.log" :::: "$TASK_LIST"
elapsed=$((SECONDS - start))

if [ "$n_todo" -gt 0 ]; then
  per_tile=$(awk -v e="$elapsed" -v n="$n_todo" -v j="$nproc" \
                 'BEGIN{printf "%.1f", e*j/n}')
  total_h=$(awk -v p="$per_tile" -v n="$N_TILES" -v j="$nproc" \
                'BEGIN{printf "%.1f", p*n/j/3600}')
  echo "[$(date +%T)] $n_todo tiles in ${elapsed}s at -j $nproc"
  echo "  ~${per_tile}s of CPU per tile -> ~${total_h}h for all $N_TILES tiles"
fi

FAILED=$(awk 'NR>1 && $7 != 0' "$WORK/parallel${SUF}.log" 2>/dev/null | wc -l)
echo "Allelome.PRO2 failures this pass: $FAILED / $n_todo"

##### ---- collect ---- #####
collect_tables final
cp "$WORK/parallel${SUF}.log" "$FINAL_DIR/" 2>/dev/null

# Count complete tables, not files that exist: the whole point of the guard is
# that those two numbers can differ.
n_final=0
for lt in "$FINAL_DIR"/*/locus_table.txt; do
  [ -e "$lt" ] || continue
  locus_table_ok "$lt" && n_final=$((n_final + 1))
done
echo "Results in $FINAL_DIR: $n_final / $N_TILES tiles have a locus table"
if [ "$n_final" -lt "$N_TILES" ] && [ "$MAX_TILES" -eq 0 ]; then
  echo "Incomplete - resubmit the same command and it will run only what is left."
fi
echo "Tile BAMs kept in $BAM_DIR - delete once you are happy with the output."
