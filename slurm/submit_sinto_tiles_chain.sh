#!/usr/bin/env bash
set -euo pipefail

# Submit spatial_sinto_tiles.slurm as a CHAIN of small, easily-scheduled jobs.
#
# Why a chain rather than one long job. Every stage of spatial_sinto_tiles
# resumes - the split, the read filter and the Allelome.PRO2 pass all skip work
# that is already complete and validated - so a job killed at the wall clock
# costs one resubmission and no recomputation. That makes a queue trade
# available: a job that sits PENDING loses that wall clock outright, whereas a
# killed one loses nothing. So where placement is slow, the chain trades job
# length for job count and comes out ahead. Where placement is immediate - as
# on serial_std - it is simply a way to run unattended past one wall clock. --dependency=afterany links them, so this is fire and
# forget - no login session, no polling (unlike submit_sc_allelome_chunks.sh,
# which has to stay attached).
#
# afterany, deliberately, NOT afterok: a link that hits the wall clock exits
# TIMEOUT, and that is the normal case here, not a failure. afterok would stop
# the chain exactly when it is most needed.
#
# One chain PER SAMPLE, not one array job of 2. An --array=1-2 link needs two
# allocations before the chain can advance; a single-task link needs one, and
# 9w and 78w then make progress independently of each other.
#
#   ./slurm/submit_sinto_tiles_chain.sh 3 64 no_Xist dup
#
# Then watch it with:
#   squeue --clusters=cm4 -u $USER -o '%.10i %.12j %.8T %.10M %.10l %.6D %R'
#
# Resources come from the #SBATCH lines at the top of spatial_sinto_tiles.slurm
# and are NOT restated here - that file is the single source of truth, so
# editing its header changes the chain too. Override only for a one-off:
#   CPUS=8 MEM=20G TIME=04:00:00 SAMPLES=1 ./slurm/submit_sinto_tiles_chain.sh 4

# 3, matching the 24h wall clock in the job header: a sample needs ~1.4 jobs of
# that length, so three links cover it with margin. Raise this if the header's
# --time is cut, and see the tables in spatial_sinto_tiles.slurm for the
# arithmetic - the two numbers have to be chosen together.
N_LINKS="${1:-3}"
TILE_UM="${2:-64}"
SNP_LABEL="${3:-no_Xist}"
FILTER_MODE="${4:-dup}"

# Empty by default, so the job's own #SBATCH lines apply. Set any of these in
# the environment to override that job for this chain only. Keep cores and
# memory a matched pair if you do - the job derives -j from both, and the
# arithmetic is documented in spatial_sinto_tiles.slurm.
CPUS="${CPUS:-}"
MEM="${MEM:-}"
TIME="${TIME:-}"
# Array task ids to chain: 1 = 9w, 2 = 78w.
SAMPLES="${SAMPLES:-1,2}"

# Seed the chain BEHIND a job that is already running, instead of starting a
# link immediately alongside it. This is not a convenience - two jobs sharing
# one $WORK is a genuine data race, because each one purges locus tables it
# judges truncated and the other may be mid-write on exactly those. Pass the
# running job's id and each per-sample chain waits on that job's MATCHING array
# task, so 9w queues behind <id>_1 and 78w behind <id>_2:
#
#   START_AFTER=209712 ./slurm/submit_sinto_tiles_chain.sh 8 64 no_Xist dup
START_AFTER="${START_AFTER:-}"
if [ -n "$START_AFTER" ] && ! [[ "$START_AFTER" =~ ^[0-9]+$ ]]; then
  echo "ERROR: START_AFTER='$START_AFTER' must be a bare numeric job id" >&2
  echo "       (the array task is appended per sample, so do not pass 209712_1)" >&2
  exit 1
fi

SCRIPT="slurm/spatial_sinto_tiles.slurm"
if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: $SCRIPT not found - run this from the Postdoc repo root." >&2
  exit 1
fi
case "$FILTER_MODE" in
  dedup|dup|none) ;;
  *) echo "ERROR: filter mode '$FILTER_MODE' - expected dedup, dup or none" >&2; exit 1 ;;
esac
if ! [[ "$N_LINKS" =~ ^[0-9]+$ ]] || [ "$N_LINKS" -lt 1 ]; then
  echo "ERROR: link count '$N_LINKS' is not a positive integer" >&2; exit 1
fi

# Only the overrides that were actually set, so an unset one leaves the job's
# own #SBATCH directive in force rather than passing an empty flag.
res=()
[ -n "$CPUS" ] && res+=(--cpus-per-task="$CPUS")
[ -n "$MEM" ]  && res+=(--mem="$MEM")
[ -n "$TIME" ] && res+=(--time="$TIME")

echo "chain: $N_LINKS links per sample"
echo "job:   $SCRIPT $TILE_UM 0 $SNP_LABEL $FILTER_MODE"
echo "tasks: $SAMPLES"
if [ ${#res[@]} -gt 0 ]; then
  echo "over:  ${res[*]}"
else
  echo "res:   from the #SBATCH header of $SCRIPT"
fi
echo

IFS=',' read -r -a task_ids <<< "$SAMPLES"
for task in "${task_ids[@]}"; do
  # An already-running array job is depended on per task, not as a whole: 9w
  # need not wait for 78w to finish.
  prev="${START_AFTER:+${START_AFTER}_${task}}"
  echo "--- array task $task ---"
  [ -n "$prev" ] && echo "  first link waits on $prev"
  for (( i=1; i<=N_LINKS; i++ )); do
    dep=()
    [ -n "$prev" ] && dep=(--dependency="afterany:$prev")
    # ${dep[@]+...}, not "${dep[@]}": under `set -u` bash before 4.4 treats an
    # empty array expansion as an unbound variable and aborts, which would kill
    # the very first link of every chain. LRZ login nodes are not guaranteed to
    # be 4.4+, so do not simplify this back.
    out=$(sbatch --parsable \
            --array="$task" \
            ${res[@]+"${res[@]}"} \
            ${dep[@]+"${dep[@]}"} \
            "$SCRIPT" "$TILE_UM" 0 "$SNP_LABEL" "$FILTER_MODE")
    # --parsable gives "<jobid>" or "<jobid>;<cluster>"; the dependency wants
    # only the id, and on a multi-cluster submit sbatch appends the cluster.
    jobid=${out%%;*}
    if [ -z "$jobid" ]; then echo "ERROR: could not parse job id from '$out'" >&2; exit 1; fi
    printf '  link %2d/%d: %s%s\n' "$i" "$N_LINKS" "$jobid" \
           "${prev:+  (after $prev)}"
    prev=$jobid
  done
  echo
done

echo "Submitted. Each link resumes where the previous one stopped; the chain"
echo "is safe to cancel at any point with: scancel -u \$USER --name=spatial_sinto_tiles"
