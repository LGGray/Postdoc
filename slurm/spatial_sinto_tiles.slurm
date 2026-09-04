#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_sinto_tiles
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
# ALTERNATIVE QUEUE, and the original one. Swap back to these three lines -
# commenting out the serial two above - if serial_std turns out to be the more
# contended of the two. serial_std takes no --qos; cm4_tiny does.
#
##SBATCH --clusters=cm4
##SBATCH --partition=cm4_tiny
##SBATCH --qos=cm4_tiny
#
# WHY serial, as of 2026-09-04. Job 209687 sat PD (Priority) on cm4_tiny with
# START_TIME N/A, i.e. the backfill scheduler could not even place it in its
# lookahead window. sprio explained why, and the explanation is not a shape
# problem at all:
#
#   JOBID  PRIORITY     AGE  FAIRSHARE  JOBSIZE  QOS
#   209687 23396878  375331   23016146     5402    0
#
# FAIRSHARE is 23.0M of the 23.4M total, and sshare put the fairshare FACTOR at
# 0.011168, so the implied weight is 23016146/0.011168 = ~2.06e9. That means the
# job was carrying 1.1% of the fairshare priority available, and a user at
# factor 1.0 would outrank it by ~88x. AGE contributes 375k against a 2.06e9
# weight - 0.02% - so WAITING CANNOT RESCUE THIS JOB, and neither can a smaller
# request: 23.4M is a fairshare verdict, not a resource one.
#
# The cause is recent usage: sshare gave EffectvUsage 0.000721 against
# NormShares 0.000172, i.e. 4.2x the normalised share, which is what the
# earlier tile runs on this very dataset consumed. That decays on the
# scheduler's half-life, so cm4_tiny recovers on its own with time - but a
# separate cluster does not have to wait for it. Confirm the weights and the
# decay before trusting any of the above:
#
#   scontrol show config | grep -i -E 'priority|backfill'
#
# The two reason codes call for opposite responses, and confusing them wastes
# days. (Resources) means the request does not fit anything free, so shrinking
# cores/memory/time is the fix - that is what the 18/40G/6h shape below is for.
# (Priority) means it fits fine but other jobs are ahead of it, and no amount of
# shrinking moves the job up that list. The levers there are different:
#
#   1. A SHORTER wall clock, which is the one shrink that still helps, because
#      it lets the backfill scheduler drop the job into a gap ahead of a large
#      reservation. A short --time with more chain links is strictly better
#      placed than a long one with fewer - the opposite of the trade that
#      applies once the job is placing immediately, which is why the 24h below
#      belongs to serial_std and would be wrong back on cm4_tiny.
#   2. A DIFFERENT QUEUE, i.e. the swap above. cm4_tiny and serial_std are
#      separate clusters with separate contention, and 18 cores / 100G / 8h is
#      already proven on serial_std by map_adult_bodymap.slurm and its siblings.
#   3. Waiting, if fairshare is the component that is low - nothing in this file
#      changes that. Check before assuming otherwise:
#        squeue --start -j <jobid>          # estimated start, the first thing to look at
#        sprio -j <jobid>                   # which priority component is short
#        sshare -U -u $USER                 # fairshare, i.e. recent usage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --array=1-2

# 18 cores / 40G / 24h. The wall clock is set by WHERE THIS RUNS, so change it
# together with the cluster block above rather than on its own.
#
# The per-tile cost is measured, not guessed. Job 207291 did 701/3738 tiles (9w)
# and 711/3215 (78w) in 2h35 at -j 30: ~6.6 min of CPU per tile, and that cost
# is the per-invocation genome-wide SNP and annotation parse, so it is fixed per
# tile and does not fall as the run proceeds. FILTER_MODE=dup scores ~2.2x the
# reads, so call it ~8.5 min per tile there. At -j 16 that is ~113 tiles/h:
#
#    6h -> ~680 tiles   9w needs 5.5 jobs
#   12h -> ~1355 tiles  9w needs 2.8 jobs
#   24h -> ~2711 tiles  9w needs 1.4 jobs, 78w needs 1.2
#
# 6h stood here while the job was queued on cm4_tiny, and for THAT queue it was
# right: pending time is the only wall clock genuinely lost when every stage
# resumes, so a short job that gets placed beats a long one that does not. On
# serial_std the job starts immediately, so that argument is gone and the trade
# inverts - a longer wall clock costs nothing and saves resubmissions. Hence
# 24h. If the cluster block is ever swapped back to cm4_tiny, cut this to 6h
# with it.
#
# 24h is also the cm4_tiny hard ceiling, so it stays valid either way. The
# serial partition's limit may well be higher, and if it is, 48h would finish
# either sample in one job. Check before raising it - sbatch rejects an
# over-limit request outright, which is harmless but confusing:
#
#   scontrol --clusters=serial show partition serial_std | tr ' ' '\n' | grep -i maxtime
#
# A long job is safe here specifically because of the resume design: a kill at
# the wall clock leaves finished locus tables in $AP2_DIR on $SCRATCH, and the
# NEXT job's startup collect_tables call rescues them, so nothing is
# recomputed. The one real hazard is leaving them there, since $SCRATCH is
# purged on an access-time policy - so always keep a follow-up link queued
# rather than ending on a killed job. The chain does that:
#
#   ./slurm/submit_sinto_tiles_chain.sh 3 64 no_Xist dup
#
# Three links of 24h rather than eight of 6h: same coverage, fewer handoffs.
# Seed it behind anything already running with START_AFTER=<jobid>; two jobs on
# one $WORK race over the locus-table purge.
#
# 40G, down from 120G, and that is a fix rather than a trim. The old note here
# recorded that "the footprint scales with -j, not with the tile size" and that
# 48G had died on SIGKILL (job 207231), and concluded the request had to be
# large. The scaling was real; the conclusion was not. The cause was GNU `sort`:
# Allelome.PRO2.sh pipes the 630 MB genome-wide SNP bed through it four times
# with no -S, and sort with no -S sizes its buffer from the NODE's physical
# memory (~512G here), not from the SLURM cgroup. Thirty workers each helped
# themselves to a share of 512G, so no request was reliably "enough" - it was a
# race against the cgroup, which 120G happened to win and 48G lost.
#
# The sort buffer is now capped explicitly (see SORT_MB below), which makes the
# per-worker footprint a chosen number: 4G for the job shell plus 2G per worker.
# 40G is that arithmetic solved for 18 cores - 4096 + 16*2048 = 36864M, so 40G
# covers the 16 workers that 18 cores allow with a little slack. The two
# numbers are a MATCHED PAIR: -j is now derived from BOTH cores and memory, so
# raising cores without raising memory leaves cores idle, and raising memory
# without raising cores buys nothing. The job prints which constraint bound it.
#
#   18 cores + 40G -> -j 16, core-limited      (this default)
#   18 cores + 32G -> -j 14, memory-limited, 2 cores idle
#   32 cores + 64G -> -j 30, core-limited      (the old shape, now placed slower)
#
# A smaller --mem is therefore safe - it lowers concurrency instead of earning a
# SIGKILL - so trimming this to get placed is a legitimate move, unlike before.
# Any change still wants a MaxRSS to point at:
#   sacct -j <jobid> -o JobID,ReqMem,MaxRSS,Elapsed,State,ExitCode
#
# The 24h cm4_tiny ceiling remains the hard cap if the wall clock is ever
# raised, but raising it is the wrong lever here: it makes placement slower for
# no gain the resume logic does not already provide.
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
#
# Arg 4 selects the read filter (see FILTER_MODE below). The duplicate-inclusive
# run wants more wall clock than the 12h above: measured on the 9w tile BAMs,
# dedup leaves 10.6 GB of the 46.8 GB split, so the MAPQ-only filter leaves
# roughly 2.2x that - ~23 GB, and ~2.2x the reads to intersect per tile.
# The dominant per-tile cost is the fixed genome-wide SNP parse, so the total
# does not scale with 2.2, but the bedtools intersect is not free either;
# ~14h becomes ~18h on that split. Give it the cm4_tiny ceiling:
#
#   sbatch --time=24:00:00 ~/Postdoc/slurm/spatial_sinto_tiles.slurm 64 0 no_Xist dup

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
# Optional arg 4: which read filter to apply to the tile BAMs before scoring.
# Mirrors the counting modes on the pysam side (spatial_tile_locus_map.slurm),
# and like SNP_LABEL it extends every output path so a second mode lands
# BESIDE the first rather than on top of it.
#
#   dedup  -F 0x400 -q 255   duplicates DROPPED, multimappers dropped. Default,
#                            and the only mode with unsuffixed paths, so every
#                            existing 64um result is byte-identical.
#   dup            -q 255    duplicates KEPT, multimappers dropped. Paths gain
#                            _dup, matching ase_pysam_dup_64um on the pysam side.
#   none   no filter at all  duplicates AND multimappers kept. Paths gain _raw.
#
#   sbatch $0 64                        # dedup, the original run
#   sbatch $0 64 0 no_Xist dup          # duplicate-inclusive, alongside it
#
# The filter suffix and the mask suffix are indistinguishable to the downstream
# R side, which builds its directory from SNP_LABEL the same way this builds
# SUF. So the dup tree is plotted by passing the filter suffix AS the label -
# no change needed there, and the PDFs come out labelled "dup":
#
#   sbatch $0 64 0 no_Xist dup                              # this job
#   sbatch ~/Postdoc/slurm/spatial_tile_map.slurm 64 9w,78w dup
#
# What "dup" is for: it is the read-level counterpart of the deduplicated run,
# so the two together separate an allelic signal from an amplification one. On
# its own it is NOT a ratio to interpret - a single over-amplified molecule can
# carry a whole tile (see the deduplication rationale below). Its use is the comparison.
FILTER_MODE="${4:-dedup}"
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
if [ "$SNP_LABEL" = "no_Xist" ]; then MASK_SUF=""; else MASK_SUF="_$SNP_LABEL"; fi

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
#
# FILTER_MODE (arg 4) chooses which half of that is applied. Everything above
# is the argument for `dedup`; `dup` keeps the duplicates deliberately, to
# quantify how much of a tile's ratio they were carrying.
FILTER_MAPQ=255
case "$FILTER_MODE" in
  dedup) FILTER_ARGS=(-F 0x400 -q "$FILTER_MAPQ"); FLT_SUF="";     FLT_SUB=bam_flt
         FILTER_DESC="duplicates, MAPQ < $FILTER_MAPQ" ;;
  dup)   FILTER_ARGS=(-q "$FILTER_MAPQ");          FLT_SUF="_dup"; FLT_SUB=bam_dup
         FILTER_DESC="MAPQ < $FILTER_MAPQ only - duplicates KEPT" ;;
  none)  FILTER_ARGS=();                           FLT_SUF="_raw"; FLT_SUB=""
         FILTER_DESC="none" ;;
  *)     echo "Unknown FILTER_MODE '$FILTER_MODE' (arg 4): expected dedup, dup or none" >&2
         exit 1 ;;
esac
# One suffix for every output path, so the mask and the filter mode are both
# free to vary without either run overwriting the other.
SUF="${MASK_SUF}${FLT_SUF}"
echo "Read filter: $FILTER_MODE  (samtools view ${FILTER_ARGS[*]:-<none>})"

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

##### ---- bound the per-worker memory, so the request can be made small ---- #####
#
# WHY 120G/32 cores used to be the only workable shape. Allelome.PRO2.sh pipes
# the 630 MB genome-wide SNP bed through GNU `sort` four times (its lines 214,
# 217, 219, 238); line 217 is the big one, since a genome-wide annotation means
# essentially every SNP survives the intersect. GNU sort with no -S sizes its
# buffer from the NODE's physical memory - ~512G on a cm4 node - and knows
# nothing about the SLURM cgroup. So 30 concurrent workers each helped
# themselves to a share of 512G and the cgroup killed the job. That, not the
# tile size, is what SIGKILLed 207231 at 48G, and it is why the footprint was
# observed to "scale with -j, not with the tile size".
#
# Capping the buffer explicitly breaks the coupling to the node: per-worker
# memory becomes a number chosen here instead of a fraction of 512G. The point
# is scheduling - a job asking 18 cores and 40G finds a slot on a shared
# partition far sooner than one asking 32 cores and 120G, and since every stage
# resumes, several small jobs that start now beat one big job that waits.
#
# Allelome.PRO2.sh is an external toolkit and is not edited. A `sort` shim early
# on PATH is how the cap reaches it. Any -S the real call already carries would
# come later on the command line and still win, so this can only lower the
# ceiling, never override an explicit choice.
# NOTE: --export=NONE at the top means `SORT_MB=512 sbatch ...` does NOT reach
# the job, for the same reason TILE_UM had to be made positional. To override at
# submit time use `sbatch --export=ALL,SORT_MB=512 ...`; otherwise edit here.
SORT_MB="${SORT_MB:-256}"
SHIM_DIR=$WORK/bin
mkdir -p "$SHIM_DIR"
# Resolve the real sort with the shim dir removed from PATH, so a re-run does
# not build a shim that execs the previous shim and loops forever.
REAL_SORT=$(PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vxF "$SHIM_DIR" | paste -sd: -) command -v sort)
if [ -z "$REAL_SORT" ]; then echo "ERROR: cannot find the real sort binary" >&2; exit 1; fi
printf '#!/bin/bash\nexec %s -S %sM "$@"\n' "$REAL_SORT" "$SORT_MB" > "$SHIM_DIR/sort"
chmod +x "$SHIM_DIR/sort"
export PATH="$SHIM_DIR:$PATH"
echo "sort shim: $SHIM_DIR/sort -> $REAL_SORT -S ${SORT_MB}M"

# nproc from BOTH cores and memory, whichever is tighter. Previously it came
# from cores alone, so a small --mem was silently a SIGKILL waiting to happen
# and the only safe way to shrink the job was not to shrink it. With the cap
# above a worker needs the sort buffer plus bedtools, awk and mpileup working
# set; 2G is that with room to spare. Raise MEM_PER_WORKER_MB, not --mem, if a
# sacct MaxRSS ever says otherwise.
MEM_PER_WORKER_MB="${MEM_PER_WORKER_MB:-2048}"
nproc=$((SLURM_CPUS_PER_TASK - 2))
if [ "$nproc" -lt 1 ]; then nproc=1; fi
# SLURM_MEM_PER_NODE is MB and is what --mem sets. --mem-per-cpu sets
# SLURM_MEM_PER_CPU instead, so honour whichever arrived.
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
  mem_mb=$SLURM_MEM_PER_NODE
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
  mem_mb=$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK))
else
  mem_mb=0
fi
if [ "$mem_mb" -gt 0 ]; then
  # 4G held back for the job shell, parallel itself and the filesystem cache.
  mem_j=$(( (mem_mb - 4096) / MEM_PER_WORKER_MB ))
  [ "$mem_j" -lt 1 ] && mem_j=1
  if [ "$mem_j" -lt "$nproc" ]; then
    echo "nproc: memory-limited to $mem_j (cores would allow $nproc)"
    echo "  ${mem_mb}M allocated, ${MEM_PER_WORKER_MB}M per worker, 4096M reserved"
    echo "  raise --mem to $(( 4096 + nproc * MEM_PER_WORKER_MB ))M to use all $nproc cores"
    nproc=$mem_j
  else
    echo "nproc: $nproc (cores), memory allows $mem_j - not limiting"
  fi
else
  echo "nproc: $nproc (cores); no memory allocation visible, not capping"
fi

##### ---- filter duplicates and multimappers ---- #####
# Mode-specific scratch directory. bam_flt stays the dedup+MAPQ output so the
# existing filtered tiles are reused untouched; bam_dup is the MAPQ-only set.
# ~10.6 GB per sample for bam_flt, so expect ~23 GB for bam_dup - the duplicates
# are ~54% of the reads and it is those bytes that come back.
if [ -n "$FLT_SUB" ]; then
  FLT_DIR=$WORK/$FLT_SUB
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
    echo "samtools view -b ${FILTER_ARGS[*]} -o '$FLT_DIR/$t' '$bam' && samtools index '$FLT_DIR/$t'"
  done >> "$WORK/filter-list.txt"
  n_flt=$(wc -l < "$WORK/filter-list.txt")
  if [ "$n_flt" -gt 0 ]; then
    echo "[$(date +%T)] filtering $n_flt tile BAMs (dropping: $FILTER_DESC)"
    parallel -j "$nproc" :::: "$WORK/filter-list.txt"
  else
    echo "[$(date +%T)] filtered tile BAMs already present"
  fi
  SCORE_DIR=$FLT_DIR
  # Worth eyeballing once. Under dedup this should be roughly a quarter of the
  # input, and the tile ratios should stop swinging between 0.1 and 1.0 once it
  # is applied. Under dup only the multimappers go, so expect most reads to
  # survive - if this drops by half, the dup flags are still being applied.
  ex=$(ls "$BAM_DIR"/*.bam | head -1); exb=$(basename "$ex")
  echo "  example $exb: $(samtools view -c "$ex") reads -> $(samtools view -c "$FLT_DIR/$exb") after filtering"
  if [ "$FILTER_MODE" = "dup" ]; then
    echo "  NOTE: duplicates are KEPT in this mode. These per-tile ratios are"
    echo "        read-level and are only meaningful next to the dedup run."
  fi
else
  SCORE_DIR=$BAM_DIR
  echo "WARNING: FILTER_MODE=none, scoring raw tile BAMs. Per-tile ratios will be"
  echo "         driven by PCR duplication and multimapping, not by allele composition."
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
# the read filter above even though the sweep counted >= SINTO_MIN_UMI chrX
# molecules, because that counter works on molecules and does not apply the
# MAPQ filter. Under FILTER_MODE=dup this should be rarer, not commoner - only
# the multimapper filter can empty a tile there.
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
