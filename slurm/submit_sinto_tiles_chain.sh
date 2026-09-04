#!/usr/bin/env bash
set -euo pipefail

# Submit spatial_sinto_tiles.slurm as a CHAIN of small, easily-scheduled jobs.
#
# Why a chain rather than one long job. Every stage of spatial_sinto_tiles
# resumes - the split, the read filter and the Allelome.PRO2 pass all skip work
# that is already complete and validated - so a job killed at the wall clock
# costs one resubmission and no recomputation. That makes a queue trade
# available: one 24h job at 32 cores and 120G can sit pending for a long time on
# a shared partition, whereas a 6h job at 16 cores and 32G starts far sooner,
# and eight of them back to back deliver more wall clock than the big job that
# is still waiting. --dependency=afterany links them, so this is fire and
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
#   ./slurm/submit_sinto_tiles_chain.sh 8 64 no_Xist dup
#
# Then watch it with:
#   squeue --clusters=cm4 -u $USER -o '%.10i %.12j %.8T %.10M %.10l %.6D %R'
#
# Overridable, with the reasoning in spatial_sinto_tiles.slurm:
#   CPUS=16 MEM=32G TIME=06:00:00 SAMPLES=1 ./slurm/submit_sinto_tiles_chain.sh 8

N_LINKS="${1:-8}"
TILE_UM="${2:-64}"
SNP_LABEL="${3:-no_Xist}"
FILTER_MODE="${4:-dup}"

# 16 cores / 32G is the shape where neither constraint is wasted: the script
# reserves 4G and budgets 2G per worker, so 32G allows exactly the 14 workers
# that 16 cores allow. Asking for more memory than that buys no parallelism,
# it only makes the job harder to place.
CPUS="${CPUS:-16}"
MEM="${MEM:-32G}"
TIME="${TIME:-06:00:00}"
# Array task ids to chain: 1 = 9w, 2 = 78w.
SAMPLES="${SAMPLES:-1,2}"

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

echo "chain: $N_LINKS links x $TIME at $CPUS cores / $MEM"
echo "job:   $SCRIPT $TILE_UM 0 $SNP_LABEL $FILTER_MODE"
echo "tasks: $SAMPLES"
echo

IFS=',' read -r -a task_ids <<< "$SAMPLES"
for task in "${task_ids[@]}"; do
  prev=""
  echo "--- array task $task ---"
  for (( i=1; i<=N_LINKS; i++ )); do
    dep=()
    [ -n "$prev" ] && dep=(--dependency="afterany:$prev")
    # ${dep[@]+...}, not "${dep[@]}": under `set -u` bash before 4.4 treats an
    # empty array expansion as an unbound variable and aborts, which would kill
    # the very first link of every chain. LRZ login nodes are not guaranteed to
    # be 4.4+, so do not simplify this back.
    out=$(sbatch --parsable \
            --array="$task" \
            --cpus-per-task="$CPUS" \
            --mem="$MEM" \
            --time="$TIME" \
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
