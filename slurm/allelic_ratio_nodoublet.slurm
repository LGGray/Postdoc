#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J allelic_ratio_nodoublet
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=24:00:00

# Re-run the allelic ratio pipeline with the scDblFinder doublet nuclei
# excluded, writing into RESULTS_ROOT=Allelic_ratio_results_nodoublet so the
# existing Allelic_ratio_results/ is untouched and the two can be diffed.
#
# No BAM is re-scored and heart_seurat_object_SCT.rds is not rebuilt: the
# doublet barcodes are dropped as each Allelome.PRO2 tree is ingested
# (drop_doublets() in 00_functions.R). Everything downstream reads
# ALLELIC_RATIOS_FILE / CEB_RATIOS_FILE under RESULTS_ROOT, so it follows.
#
# 24h is the serial_std walltime cap. 4 cpus / 100G is the allocation the other
# jobs over this object (lox_calling, allelic_ratio_cutoff_sweep) run in; memory
# is dominated by readRDS on the 5 GB object, not by the stages.
#
# With the FPR simulation off by default (see FPR_N_REPS) the long pole is gone
# and the chain should finish well inside the cap; STAGES below still lets it be
# split if some stage surprises us.

module load slurm_setup

# --export=NONE means STAGES=... on the sbatch line does not reach the job.
# Positional arguments instead:
#
#   sbatch slurm/allelic_ratio_nodoublet.slurm                    # full chain
#   sbatch slurm/allelic_ratio_nodoublet.slurm "12 02 04"         # ingest + core
#   sbatch slurm/allelic_ratio_nodoublet.slurm "03 05 06 07 08 11"  # the rest
#   sbatch slurm/allelic_ratio_nodoublet.slurm "12 02" 200        # cheap FPR pass
#
# $1  stages to run, in order. Default is the full dependency-ordered chain:
#       12 ingest the EXISTING table -> whole_chr_allelic_ratios.txt
#          12, not 10: the read-level Allelome.PRO2 tree that 10 walks is no
#          longer on disk. 12 filters the table 10 previously produced, and
#          also corrects allelic_ratio from ar_a1 to ar_dom - see the header of
#          12_filter_existing_table.R. Switch back to 10 if that tree is ever
#          rebuilt with slurm/scAllelome.slurm; the per-cell BAMs are intact.
#       02 whole chrX           -> whole_chr_cell_metadata.txt  (needs 12)
#       04 core escape block    -> core_escape_block_*.txt
#       03 per gene             (needs 02)
#       05 LOX sensitivity      (needs 02 and 04)
#       06 cutoff sweep         (needs 12)
#       07 core escape sweep    (needs 04)
#       08 LOX calling          (needs 04)
#       11 autosomal control    (needs 12)
# $2  FPR_N_REPS for 02's dispersion false-positive simulation. DEFAULT 0,
#     i.e. skipped -- that simulation exists to quantify the n=1-per-condition
#     caveat, and it already did: the dispersion LRTs run at 8-81% false
#     positive rate, so they are reported descriptively and not as tests.
#     Re-running ~48,000 glmmTMB fits on the no-doublet cells would re-derive a
#     conclusion that does not depend on the doublets. The existing
#     Allelic_ratio_results/cutoff_30/whole_chr_*_dispersion_FPR.txt stays the
#     reference for that caveat. Pass 200 or 1000 here to run it anyway.
STAGES="${1:-12 02 04 03 05 06 07 08 11}"
FPR_N_REPS="${2:-0}"

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, not the RNAseq env most of slurm/ uses: these are Seurat scripts.
conda activate seurat_env

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
cd $BASE/OCM

# The scDblFinder per-nucleus calls, written by OCM_heart/scDblFinder_rates.R.
export DOUBLET_FILE=Allelic_ratio_results/scDblFinder_per_cell.txt
export RESULTS_ROOT=Allelic_ratio_results_nodoublet
export FPR_N_REPS

if [ ! -f "$DOUBLET_FILE" ]; then
  echo "Missing $DOUBLET_FILE - run slurm/scDblFinder_rates.slurm first." >&2
  exit 1
fi

# Fail before walking any tree if the classes are not where 00_functions.R
# expects them, rather than an hour into stage 10.
head -1 "$DOUBLET_FILE" | grep -q "class" || {
  echo "$DOUBLET_FILE has no 'class' column" >&2; exit 1; }

echo "DOUBLET_FILE: $DOUBLET_FILE"
echo "RESULTS_ROOT: $RESULTS_ROOT"
echo "STAGES:       $STAGES"
echo "FPR_N_REPS:   $FPR_N_REPS"

# A case rather than an associative array: `declare -A` needs bash >= 4.0, and
# this repo already carries a bash < 4.4 workaround in
# submit_sinto_tiles_chain.sh, so the shell here is not assumed to be new.
stage_script() {
  case "$1" in
    12) echo 12_filter_existing_table.R ;;
    10) echo 10_build_ratio_table.R ;;
    02) echo 02_whole_chrX.R ;;
    03) echo 03_all_genes.R ;;
    04) echo 04_core_escape.R ;;
    05) echo 05_lox_sensitivity.R ;;
    06) echo 06_cutoff_sweep.R ;;
    07) echo 07_core_escape_cutoff_sweep.R ;;
    08) echo 08_lox_calling.R ;;
    11) echo 11_autosomal_control.R ;;
    *)  echo "" ;;
  esac
}

# Validate the whole list before running anything: a typo'd stage should not
# surface three hours in, after 10 and 02 have already been paid for.
for s in $STAGES; do
  if [ -z "$(stage_script "$s")" ]; then
    echo "Unknown stage '$s'. Valid: 12 10 02 03 04 05 06 07 08 11" >&2
    exit 1
  fi
done

for s in $STAGES; do
  f="$(stage_script "$s")"
  echo ""
  echo "=== stage $s : $f ==============================================="
  date
  # Deliberately not `set -e`: a late stage failing should not throw away the
  # earlier ones, and 03/05 can legitimately fail on a tree that lacks the
  # per-gene or core-escape annotation. Each stage's status is reported instead.
  if Rscript "$BASE/Postdoc/OCM_heart/allelic_ratio/$f"; then
    echo "--- stage $s OK"
  else
    echo "--- stage $s FAILED (exit $?) - continuing" >&2
  fi
done

echo ""
echo "Done. Results in $PWD/$RESULTS_ROOT"
if [ "$FPR_N_REPS" -le 0 ]; then
  echo ""
  echo "NOTE: the dispersion LRT tables under $RESULTS_ROOT were written, but"
  echo "      their FPR calibration was skipped (FPR_N_REPS=0). With n=1 animal"
  echo "      per condition those p-values are descriptive, not tests -- see"
  echo "      Allelic_ratio_results/cutoff_*/whole_chr_*_dispersion_FPR.txt."
fi
