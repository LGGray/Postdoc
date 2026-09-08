#!/bin/bash
# Submit the overnight Allelome.PRO2 run and chain the plotting job behind it.
#
#   bash slurm/submit_allelome_chain.sh              # submit both
#   bash slurm/submit_allelome_chain.sh 5476689      # chain onto an EXISTING run job
#
# afterany, deliberately: the run job's per-nucleus stage may not finish inside
# 24h, and partial results are still worth plotting. The per-celltype stage runs
# first precisely so something usable exists either way.
#
# THE JOB-ID STRIP IS REQUIRED. On this multi-cluster setup `sbatch --parsable`
# returns "<jobid>;<cluster>", e.g. "5476689;serial", not a bare number. Passing
# that straight into --dependency gives afterany:5476689;serial and sbatch
# rejects it with "Job dependency problem" - which looks like a scheduler
# complaint about the dependency itself rather than a parsing mistake.
set -euo pipefail
cd "$(dirname "$0")/.."

if [ $# -ge 1 ]; then
  RUN_ID="${1%%;*}"
  echo "chaining onto existing run job: $RUN_ID"
else
  RUN_RAW=$(sbatch --parsable --clusters=serial slurm/multiome_allelome_run.slurm)
  RUN_ID="${RUN_RAW%%;*}"
  echo "run job:  $RUN_ID   (sbatch returned '$RUN_RAW')"
fi

if ! [[ "$RUN_ID" =~ ^[0-9]+$ ]]; then
  echo "could not parse a numeric job id from '${RUN_RAW:-$1}'" >&2
  exit 1
fi

PLOT_RAW=$(sbatch --parsable --clusters=serial \
             --dependency=afterany:"$RUN_ID" slurm/multiome_allelome_plots.slurm)
echo "plot job: ${PLOT_RAW%%;*}   (afterany:$RUN_ID)"
echo
echo "Watch with:  squeue --clusters=serial -u \$USER"
echo "If the run times out, resubmit it - every stage skips work already done -"
echo "then rerun the plot job with REBUILD=1 to rescan for the new nuclei."
