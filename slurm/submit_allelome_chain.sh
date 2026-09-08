#!/bin/bash
# Submit the overnight Allelome.PRO2 run and chain the plotting job behind it.
#
#   bash slurm/submit_allelome_chain.sh
#
# afterany, deliberately: the run job's per-nucleus stage may not finish inside
# 24h, and partial results are still worth plotting. The per-celltype stage runs
# first precisely so something usable exists either way.
set -euo pipefail
cd "$(dirname "$0")/.."

RUN=$(sbatch --parsable slurm/multiome_allelome_run.slurm)
echo "run job:  $RUN"
PLOT=$(sbatch --parsable --dependency=afterany:"$RUN" slurm/multiome_allelome_plots.slurm)
echo "plot job: $PLOT  (afterany:$RUN)"
echo
echo "Watch with:  squeue --clusters=serial -u \$USER"
echo "If the run times out, resubmit it - every stage skips work already done -"
echo "then rerun the plot job with REBUILD=1 to rescan for the new nuclei."
