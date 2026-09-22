#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_atac_plots
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=6:00:00

# Consolidate the allelic ATAC tree and plot it, including the RNA/ATAC
# concordance. Step 7 of multiome/ANALYSIS_PLAN.md.
#
#   sbatch slurm/multiome_atac_plots.slurm
#
# ---------------------------------------------------------------------------
# READS DSS, NOT SCRATCH - AND THAT IS THE WHOLE REASON THIS NEEDS A COPY FIRST
# ---------------------------------------------------------------------------
# multiome_allelome_atac.slurm writes every stage to $SCRATCH and has NO
# copy-back stage, so its output only reaches
# adult_aged_multiome/multiome_allelome_atac/ when it is rsynced there by hand:
#
#   rsync -a $SCRATCH/multiome_allelome_atac/ \
#     $BASE/adult_aged_multiome/multiome_allelome_atac/
#
# 08_atac_allelome.R therefore defaults to the DSS copy rather than to scratch,
# which is the opposite of 05_allelome_plots.R. To run against a tree still on
# scratch, uncomment the export below - it has to be set HERE and not in the
# submitting shell, because --export=NONE means nothing from that shell arrives.
#
# export ATAC_TREE=$SCRATCH/multiome_allelome_atac
#
# ---------------------------------------------------------------------------
# WHY 6h RATHER THAN THE 3h multiome_allelome_plots.slurm ASKS FOR
# ---------------------------------------------------------------------------
# The cost here is metadata, not compute: the consolidation stats and opens one
# locus_table.txt per output directory, and the per-nucleus tree is 1,422 (9w)
# plus 3,125 (78w) directories on a parallel filesystem where that is the slow
# operation. The cached TSVs make every later run trivial - the first one is
# what needs the headroom. REBUILD=1 forces a rescan past the cache.
#
# 4 cpus is not a throttle: the R is single-threaded and the extra cores are
# there for the BLAS and for read_tsv. Do not raise it expecting a speedup, and
# note serial_std caps cpus-per-task at 16 and sbatch REJECTS an over-request
# rather than trimming it.
#
# Depends on 05_allelome_plots.R having written figures_allelome/, which is
# where the RNA side of both concordance figures comes from. If that is absent
# the ATAC figures still get written and the concordance is skipped with a
# message, so this is a soft dependency and afterany is fine:
#
#   PLOTS=$(sbatch --parse slurm/multiome_allelome_plots.slurm)
#   sbatch --dependency=afterany:${PLOTS%%;*} slurm/multiome_atac_plots.slurm
# ---------------------------------------------------------------------------

module load slurm_setup
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate seurat_env
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc || exit 1
Rscript multiome/08_atac_allelome.R
