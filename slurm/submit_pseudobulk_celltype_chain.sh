#!/usr/bin/env bash
set -euo pipefail

# Submit the whole per-cell-type pseudobulk allelic-ratio pipeline as ONE
# dependency chain, so it is fire and forget - no login session, no polling
# between stages. See OCM_heart/PSEUDOBULK_CELLTYPE.md for what each stage
# does and why.
#
#   ./slurm/submit_pseudobulk_celltype_chain.sh
#   ./slurm/submit_pseudobulk_celltype_chain.sh merge Allelic_ratio_results_nodoublet
#
#   $1  mode          merge | sinto   (default merge; see the BAM job header)
#   $2  RESULTS_ROOT  default Allelic_ratio_results
#   $3  TERMINAL_MB   default 5
#
# Four stages, in order:
#   ids       barcode -> cell type maps, doublets dropped   (minutes)
#   bams      one pseudobulk BAM per (cohort, cell type)    (array 1-4, hours)
#   allelome  Allelome.PRO2 over them                       (array 1-4, the long one)
#   analysis  tables, matched contrast, figures             (minutes)
#
# afterok, NOT afterany - the opposite choice to
# submit_sinto_tiles_chain.sh, and for the opposite reason. Nothing here
# resumes: the stages are strictly sequential transformations, so a failed or
# timed-out BAM build means the Allelome.PRO2 stage would score a partial cell
# type and produce a number that looks fine and is wrong. Stop instead.
#
# Re-run one stage without redoing the rest - this is the normal way to widen
# the terminal window after reading the coverage table:
#   STAGES=analysis ./slurm/submit_pseudobulk_celltype_chain.sh merge Allelic_ratio_results 10
#
# Watch it with:
#   squeue -M cm4 -u $USER -o '%.10i %.28j %.8T %.10M %.10l %.6D %R'

MODE="${1:-merge}"
RESULTS_ROOT="${2:-Allelic_ratio_results}"
TERMINAL_MB="${3:-5}"
N_PERM="${N_PERM:-2000}"
STAGES="${STAGES:-ids,bams,allelome,analysis}"

case "$MODE" in
  merge|sinto) ;;
  *) echo "ERROR: mode '$MODE' - expected merge or sinto" >&2; exit 1 ;;
esac
if ! [[ "$TERMINAL_MB" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "ERROR: TERMINAL_MB '$TERMINAL_MB' is not a number" >&2; exit 1
fi

R_JOB="slurm/pseudobulk_celltype_R.slurm"
BAM_JOB="slurm/pseudobulk_celltype_bams.slurm"
ALL_JOB="slurm/pseudobulk_celltype_allelome.slurm"
for f in "$R_JOB" "$BAM_JOB" "$ALL_JOB"; do
  [[ -f "$f" ]] || { echo "ERROR: $f not found - run this from the Postdoc repo root." >&2; exit 1; }
done

# Every stage name given must be a real one, or a typo silently submits
# nothing and looks like success.
IFS=',' read -r -a want <<< "$STAGES"
for st in "${want[@]}"; do
  case "$st" in
    ids|bams|allelome|analysis) ;;
    *) echo "ERROR: unknown stage '$st' in STAGES='$STAGES'" >&2
       echo "       expected some comma-separated subset of: ids,bams,allelome,analysis" >&2
       exit 1 ;;
  esac
done
wants() { local s; for s in "${want[@]}"; do [[ "$s" == "$1" ]] && return 0; done; return 1; }

echo "mode:         $MODE"
echo "RESULTS_ROOT: $RESULTS_ROOT"
echo "TERMINAL_MB:  $TERMINAL_MB   N_PERM: $N_PERM"
echo "stages:       $STAGES"
echo

# --parsable gives "<jobid>" or "<jobid>;<cluster>"; a dependency wants only
# the id, and these all submit to cm4, so sbatch appends the cluster.
submit() {
  local desc="$1"; shift
  local dep=()
  [ -n "$prev" ] && dep=(--dependency="afterok:$prev")
  # ${dep[@]+...}, not "${dep[@]}": under `set -u` bash before 4.4 treats an
  # empty array expansion as unbound and aborts, which would kill the first
  # stage of every chain. LRZ login nodes are not guaranteed to be 4.4+.
  local out
  out=$(sbatch --parsable ${dep[@]+"${dep[@]}"} "$@")
  local jobid=${out%%;*}
  [ -n "$jobid" ] || { echo "ERROR: could not parse a job id from '$out'" >&2; exit 1; }
  printf '  %-9s %s%s\n' "$desc" "$jobid" "${prev:+  (afterok $prev)}"
  prev="$jobid"
}

# A dependency cannot cross LRZ clusters, which is why pseudobulk_celltype_R
# is on cm4 like the other two rather than on serial_std.
prev=""
wants ids      && submit "ids"      "$R_JOB"   ids "$RESULTS_ROOT"
wants bams     && submit "bams"     "$BAM_JOB" "$MODE"
wants allelome && submit "allelome" "$ALL_JOB"
wants analysis && submit "analysis" "$R_JOB"   analysis "$RESULTS_ROOT" "$TERMINAL_MB" "$N_PERM"

if [ -z "$prev" ]; then
  echo "  nothing submitted - STAGES='$STAGES' selected no stage" >&2
  exit 1
fi

echo
echo "Submitted. The last stage writes:"
echo "  OCM/$RESULTS_ROOT/pseudobulk_celltype/pseudobulk_celltype_chromosome_ends.pdf"
echo
echo "Cancel the lot with: scancel -M cm4 -u \$USER --name=pseudobulk_celltype_R,pseudobulk_celltype_bams,pseudobulk_celltype_allelome"
